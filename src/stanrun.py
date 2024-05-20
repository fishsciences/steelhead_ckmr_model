import cmdstanpy as stan
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import chain

def transitive_closure(data):
    ids = data['SampleID_GIQ'].unique()
    dups = data.loc[~data['DuplicateOf'].isna(), ['SampleID_GIQ', 'DuplicateOf']]
    dups = dict(zip(dups['SampleID_GIQ'].values, dups['DuplicateOf'].values))
    eqs = dict((x, i) for (i, x) in enumerate(ids))
    changed = len(eqs)
    while changed > 0:
        changed = 0
        for (l, r) in dups.items():
            if eqs[l] > eqs[r]:
                eqs[l] = eqs[r]
                changed += 1
    # labels are now unique per equivalence class, 
    # but they are not consecutive
    res = defaultdict(set)
    for (x, i) in eqs.items():
        res[i].add(x)
    res = sorted(((i, x) for (i, x) in res.items()), key=(lambda x: x[0]))
    res = enumerate(x for (_, x) in res)
    res = chain.from_iterable(((x, i) for x in xs) for (i, xs) in res)
    res = dict(res)
    res = data.assign(FishID=data['SampleID_GIQ'].apply(lambda x: res[x]))
    return res


def age_cleanup(data):
    age = data['Scale_age'].apply(lambda x: str(x).strip('+'))
    age = age.apply(lambda x: float(x) if x.isdigit() else np.nan)
    res = data.assign(Age=age)
    return res

def impute_age(meta, pairs):
    df = pairs.groupby('D1_indiv_parent')['BirthYear_offspring'].min().reset_index()
    df = df.rename(columns={
        'D1_indiv_parent': 'SampleID_GIQ',
        'BirthYear_offspring': 'BirthYear_UB'
    })
    df['BirthYear_UB'] -= 1
    res = meta.merge(df, how='left', left_on='SampleID_GIQ', right_on='SampleID_GIQ')
    res['Age_LB'] = res['CollectionYear'] - res['BirthYear_UB']
    df = res.loc[:, ['Age', 'Age_LB']].bfill(axis=1)
    res['Age'] = df['Age']
    res = res.drop(columns=['Age_LB', 'BirthYear_UB'])
    res.loc[res['lifestage'] == 'juvenile', 'Age'] = 0.0
    return res    

def verify_age(data):
    df = data.loc[:, ['SampleID_GIQ', 'DuplicateOf', 'CollectionYear', 'Age']]
    df = df.merge(df, how='inner', left_on='DuplicateOf', right_on='SampleID_GIQ')
    df['AgeConflict'] =  (~df['Age_x'].isna()) 
    df['AgeConflict'] &= (~df['Age_y'].isna()) 
    df['AgeConflict'] &= (~df['CollectionYear_x'].isna()) 
    df['AgeConflict'] &= (~df['CollectionYear_y'].isna()) 
    df['AgeConflict'] &= (df['Age_x'] - df['Age_y']) != (df['CollectionYear_x'] - df['CollectionYear_y'])
    df = df.loc[:, ['SampleID_GIQ_x', 'AgeConflict']]
    df = df.rename(columns={'SampleID_GIQ_x': 'SampleID_GIQ'})
    df = data.merge(df, how='left', left_on='SampleID_GIQ', right_on='SampleID_GIQ')
    df['AgeConflict'] = df['AgeConflict'].fillna(value=False)
    return df


def relabel_pairs(meta, pairs):
    lbl = meta.loc[:, ['SampleID_GIQ', 'FishID']]
    df = pairs.merge(lbl, how='inner', left_on='D2_indiv_offspring', right_on='SampleID_GIQ')
    df = df.rename(columns={'FishID': 'offspring'})
    df = df.merge(lbl, how='inner', left_on='D1_indiv_parent', right_on='SampleID_GIQ')
    df = df.rename(columns={'FishID': 'parent'})
    df = df.loc[:, ['parent', 'offspring', 'BirthYear_offspring']]
    df = df.rename(columns={'BirthYear_offspring': 'birthyear'})
    return df


def load_data(meta = '../data/meta_clean.csv', pairs = '../data/pairs_clean.csv'):
    meta = pd.read_csv(meta)
    pairs = pd.read_csv(pairs)
    return (meta, pairs)

def mk_stan(meta, pairs):
    Y = meta['CollectionYear'].nunique()
    M = meta['lifestage'].nunique()
    A = (1 + meta['Age'] + meta['CollectionYear'].max() - meta['CollectionYear']).max()
    A = int(A)
    Obs = len(meta)
    P = len(pairs)
    F = meta['FishID'].nunique()
    res = {'Y': Y, 'M': M, 'A': A, 'Obs': Obs, 'P': P, 'F': F}
    sex = meta.groupby('FishID')['Sex_final'].first().reset_index()
    sex = sex.sort_values(by='FishID')['Sex_final']
    sex = sex.apply(lambda x: mk_stan.sex_enc[x])
    res['fish_sex'] = sex
    res['obs_fish'] = (meta['FishID'] + 1).values
    res['obs_method'] = meta['lifestage'].apply(lambda x: mk_stan.meth_enc[x]).values
    res['obs_year'] = (meta['CollectionYear'] - meta['CollectionYear'].min() + 1).values
    res['obs_age'] = (meta['Age'] + 1).astype(np.int64).values
    res['parent'] = (pairs['parent'] + 1).values
    res['offspring'] = (pairs['offspring'] + 1).values
    res['app'] = np.zeros((M, 2, A), dtype=np.int64)
    res['app'][0, :, 0]  = 1
    res['app'][1, :, 1:] = 1
    return res

mk_stan.sex_enc  = {'male': 1, 'female': 2}
mk_stan.meth_enc = {'juvenile': 1, 'adult': 2}


def mk_stan_single(meta, pairs, mu=0.5, kappa=2.0):
    years = meta['CollectionYear'].sort_values().unique()
    res = dict()
    for y in years:
        ad = (meta['CollectionYear'] == y) & (meta['lifestage'] == 'adult')
        fs = ad & (meta['Sex_final'] == 'female')
        ms = ad & (meta['Sex_final'] == 'male')
        js = (meta['CollectionYear'] == y) & (meta['lifestage'] == 'juvenile')
        fs = meta.loc[fs] # year's females
        ms = meta.loc[ms] # year's males
        js = meta.loc[js] # year's juveniles
        ss = pd.concat([fs, ms]) # year's spawners
        sex = meta.loc[ad, ['FishID', 'Sex_final']]
        sex = sex.groupby('FishID')['Sex_final'].first().reset_index()
        adults = set(fs['FishID'].unique()).union(ms['FishID'].unique())
        juveniles = set(js['FishID'].unique())
        ps = pairs.loc[pairs['parent'].isin(adults) & pairs['offspring'].isin(juveniles)]
        ps = ps.merge(sex, how='inner', left_on='parent', right_on='FishID')
        ps = ps.groupby(['offspring', 'Sex_final'])['parent'].nunique().reset_index()
        ps = ps.loc[:, ['offspring', 'Sex_final', 'parent']]
        ps = ps.rename(columns={'Sex_final': 'sex', 'parent': 'count'})
        ps['count'] = ps['count'].apply(lambda x: min(1, x))
        ps['sex'] = ps['sex'].apply(lambda x: 1 if x == 'female' else 2)
        ps['count'] = ps['count'] * ps['sex']
        ps = ps.drop(columns=['sex'])
        ps = ps.groupby('offspring')['count'].sum().reset_index()
        ps = ps.groupby('count')['offspring'].nunique().reset_index()
        ps = ps.rename(columns={'count': 'class', 'offspring': 'count'})
        ps = ps.sort_values(by='class')
        M = ms['FishID'].nunique()
        F = fs['FishID'].nunique()
        J = js['FishID'].nunique()
        by_parent = np.zeros(4, dtype=np.int64)
        for (i, r) in ps.iterrows():
            by_parent[r['class']] = r['count']
        by_parent[0] = J - by_parent[1:].sum()
        stan_data = {'F': F, 'M': M, 'mu': mu, 'kappa': kappa, 'by_parents': by_parent}
        res[y] = stan_data
    return res

def plot_fom(smp, ci=0.95):
    alpha = 0.5 * (1 - ci)
    qs = [alpha, 1 - alpha]
    foms = [c for c in smp.columns if c.startswith('log_fom')]
    smp = smp.loc[:, foms]
    if len(foms) % 2 != 0:
        raise RuntimeError('Invalid posterior sample')
    A = len(foms) // 2
    fig, ax = plt.subplots(nrows=2, ncols=1)
    for s in range(2):
        cs = [f'log_fom[{a+1},{s+1}]' for a in range(A)]
        cs = [np.exp(smp[c].values) for c in cs] 
        ax[s].violinplot(cs, showmeans=True, showextrema=False, quantiles=len(cs)*[qs])
        ax[s].set_xlabel('Age group')
        title = f'Probability of surviving another year for {"males" if s == 0 else "females"}'
        ax[s].set_title(title)
    return (fig, ax)
