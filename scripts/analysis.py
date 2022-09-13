# import taskdb2
import pandas as pd
import nscomplex
import networkx as nx
import snappy
import json
import regina
from nscomplex import reconstruct_faces
from sage.all import (PolynomialRing, QQ, ZZ, matrix, vector,
                      gcd, lcm, SimplicialComplex, Simplex,
                      cyclotomic_polynomial, PowerSeriesRing,
                      divisors, moebius)

import collections

def dim_faces(max_faces):
    return tuple(sorted([C['dim'] for C in max_faces]))

def examine_nonnormal_iso(dataframe):
    df = dataframe[dataframe.nonnormal_iso==1].copy()
    df.max_faces = df.max_faces.apply(eval)
    return df.max_faces.apply(dim_faces).value_counts()

def examine_nonnormal_iso_carefully(dataframe):
    df = dataframe[dataframe.nonnormal_iso==1].copy()
    for i, row in df.iterrows():
        CS = nscomplex.EssentialSurfaces(row['tri_used'], row['LW_euler_bound'])
        faces = CS.essential_faces_of_normal_polytope()
        for C in faces.maximal:
            if C.dim > 0:
                assert C.dim == 1
                assert CS.LG.has_edge(*C.vertex_surfaces)
                print(row['name'])


def examine_LW(row):
    new_info = dict()
    vertex_surfaces = eval(row['vertex_surfaces'])
    vertex_genera = eval(row['vertex_genera'])
    max_faces = eval(row['max_faces'])
    all_faces = eval(row['all_faces'])

    one_skel = nx.Graph()
    one_skel.add_nodes_from(vertex_surfaces.keys())
    one_skel.add_edges_from([F['verts'] for F in all_faces if F['dim']==1])

    new_info['LW_comps'] = nx.number_connected_components(one_skel)

    # Now divide the maximal faces up

    one_skel_comps = [frozenset(verts) for verts in
                      nx.connected_components(one_skel)]
    components_to_max_faces = {verts:[] for verts in one_skel_comps}
    components_to_all_faces = {verts:[] for verts in one_skel_comps}

    for face in max_faces:
        verts = set(face['verts'])
        matches = [comp for comp in components_to_max_faces if verts.issubset(comp)]
        assert len(matches) == 1
        components_to_max_faces[matches[0]].append(face)


    for face in all_faces:
        verts = set(face['verts'])
        matches = [comp for comp in components_to_all_faces if verts.issubset(comp)]
        assert len(matches) == 1
        components_to_all_faces[matches[0]].append(face)

    new_info['LW_verts_comp'] = repr(sorted([len(comp) for comp in one_skel_comps]))

    eulers = []
    for components, faces in components_to_all_faces.items():
        euler = 0
        for face in faces:
            euler += (-1)**face['dim']
        eulers.append(euler)

    eulers = sorted(eulers)
    assert sum(eulers) == row['LW_euler']
    new_info['LW_euler_comps'] = repr(eulers)


    new_info['LW_inhom_comp'] = 0
    comp_dims = []
    for component, faces in components_to_max_faces.items():
        dims = {face['dim'] for face in faces}
        if len(dims) > 1:
            new_info['LW_inhom_comp'] = 1
        comp_dims.append(max(dims))

    new_info['LW_dim_comps'] = repr(sorted(comp_dims))

    new_info['LW_max_vert_face'] = max([len(face['verts']) for face in max_faces])

    return new_info


def add_LW_info():
    db = taskdb2.ExampleDatabase('surface_counts_cusped_new')
    df = db.dataframe()
    df = df[df.nonnormal_iso==0]
    for i, row in df.iterrows():
        new_info = examine_LW(row)
        for key, val in new_info.items():
            df.at[i, key] = val

    updated_cols = ['LW_comps', 'LW_dim_comps', 'LW_euler_comps',
                    'LW_inhom_comp', 'LW_max_vert_face', 'LW_verts_comp']
    return df, updated_cols

def draw_fancy(row):
    """
    6768   K15n27229  ...  [4,8,4,9,8,16,12,24,12,39,20,36]
    """
    from sage.all import matrix, vector, QQ
    import numpy as np
    #row = dataframe.loc[14011]
    vertex_surfaces = eval(row['vertex_surfaces'])
    vertex_genera = eval(row['vertex_genera'])
    all_faces = eval(row['all_faces'])

    one_skel = nx.Graph()
    one_skel.add_nodes_from(vertex_surfaces.keys())
    one_skel.add_edges_from([F['verts'] for F in all_faces if F['dim']==1])

    surf_vecs = []
    surf_to_index = dict()
    for i, key in enumerate(vertex_surfaces.keys()):
        surf_to_index[key] = i
        Fvec = vector(QQ, vertex_surfaces[key])
        g = vertex_genera[key]
        surf_vecs.append(Fvec/(2*g - 2))


    surf_mat = matrix(surf_vecs)
    surf_mat = surf_mat.denominator() * surf_mat
    good_basis = [r for r in surf_mat.LLL().rows() if r != 0]
    B = matrix(good_basis).transpose()
    #B = surf_mat.transpose()

    proj_surf = np.array([B.solve_right(r) for r in surf_mat.rows()],
                            dtype=np.float64)

    # center data at 0, do PCA
    proj_surf += -proj_surf.mean(axis=0)
    covar = np.dot(proj_surf.T, proj_surf)
    eigen_vals, eigen_vecs = np.linalg.eig(covar)
    surf_pca = np.dot(proj_surf, eigen_vecs)[:, :3]

    two_faces = []
    three_faces = []
    for face in all_faces:
        verts = face['verts']
        if face['dim'] == 2:
            C = one_skel.subgraph(verts)
            cycle = [u for u, v in nx.find_cycle(C)]
            assert set(cycle) == set(verts)
            two_faces.append([surf_to_index[name] for name in cycle])
        if face['dim'] == 3:
            three_faces.append([surf_to_index[name] for name in verts])


    return two_faces, three_faces, surf_pca, surf_to_index


def compare_tri_without_isos():
    db = taskdb2.ExampleDatabase('surface_counts_knots')
    df = db.dataframe()
    with_iso= pd.read_csv('data_LW/knots_nonnormal_iso.csv')
    without_iso = df[df.name.isin(set(with_iso.name))].reset_index()
    return with_iso, without_iso


# K15n142389

def save_deep_genus():
    """
    sage: dg = save_deep_genus()
    sage:
    50      9
    100    21
    200    43
    """
    db = taskdb2.ExampleDatabase('surface_counts_deep')
    df = db.dataframe()
    def eval_if_possible(data):
        if data is not None:
            data = eval(data)
        return data

    for col in ['by_genus', 'by_genus_50', 'by_genus_100', 'by_genus_200']:
        df[col] = df[col].apply(eval_if_possible)

    dg = df[df.by_genus_100.notnull() & df.by_genus_200.notnull()]
    assert all(dg.by_genus_100.apply(lambda x:x[-1]) == dg.by_genus_200.apply(lambda x:x[0]))

    def best(row):
        if row.by_genus_200 is not None and row.by_genus_100 is not None:
            return row.by_genus_100 + row.by_genus_200[1:]
        elif row.by_genus_100 is not None:
            return row.by_genus_100
        elif row.by_genus_50 is not None:
            return row.by_genus_50
        else:
            return row.by_genus

    def nonzero(L):
        return len(L) - L.count(0)

    df['by_genus'] = df.apply(best, axis=1)

    df = df[df.by_genus.apply(nonzero) > 3]

    cols = ['name',
            'volume',
            'tets',
            'nonnormal_iso',
            'LW_dim',
            'LW_euler',
            'LW_comps',
            'LW_dim_comps',
            'LW_euler_comps',
            'LW_inhom_comp',
            'LW_max_vert_face',
            'LW_verts_comp',
            'num_max_faces',
            'tri_used',
            'vertex_surfaces',
            'vertex_genera',
            'max_faces',
            'all_faces',
            'LW_euler_bound',
            'gen_func',
            'by_genus']

    df.by_genus = df.by_genus.apply(lambda x:repr(x).replace(' ', ''))
    df = df[cols]

    stats = df.by_genus.apply(lambda x:len(eval(x))).value_counts().sort_index()
    print(stats)
    return df

def overlap():
    """
    Overlap between the two samples is just 216 manifolds: 209 small
    and 7 very large.  The very large ones are:

    ['t12198', 't12200', 't12756', 't12757', 'o9_43566', 'o9_43608', 'o9_43609']
    """

    dc = taskdb2.ExampleDatabase('surface_counts_cusped_new').dataframe()
    dk = taskdb2.ExampleDatabase('surface_counts_knots_new').dataframe()

    def lookup_cusped(name):
        M = snappy.Manifold(name)
        return not snappy.OrientableCuspedCensus.identify(M) is False

    def lookup_knots(name):
        M = snappy.Manifold(name)
        if M.homology().elementary_divisors() != [0]:
            return False
        return not snappy.HTLinkExteriors.identify(M) is False

    #assert sum(dc.name.apply(lookup_knots)) == sum(dk.name.apply(lookup_cusped))

    is_dup = dc.name.apply(lookup_knots)
    print(len(dc), len(dk), len(dc) + len(dk) - len(dc[is_dup]))


def save_very_large():
    dc = taskdb2.ExampleDatabase('surface_counts_cusped_new').dataframe()
    dk = taskdb2.ExampleDatabase('surface_counts_knots_new').dataframe()
    dupes = ['t12198', 't12200', 't12756', 't12757', 'o9_43566', 'o9_43608', 'o9_43609']
    dcvl = dc[(dc.LW_dim>0)&(~dc.name.isin(dupes))]
    dkvl = dk[(dk.nonnormal_iso==0)&(dk.LW_dim>0)]
    dans = pd.concat([dcvl, dkvl], axis=0)
    dans.to_csv('./data_LW/very_large_combined.csv')

R = PolynomialRing(ZZ, 'x')
x = R.gen()
F = R.fraction_field()


def load_combined():
    df = pd.read_csv('./data_LW/very_large_combined.csv')
    df.gen_func = df.gen_func.apply(F)
    df['gen_func_denom'] = df.gen_func.apply(lambda g:g.denominator().factor())

    def x_minus_one_factor(factors):
        for p, e in factors:
            if p == x - 1:
                return e

    def gen_func_deg(gen_func):
        n = gen_func.numerator().degree()
        d = gen_func.denominator().degree()
        assert n == d
        return n

    def period(gen_func):
        factors = gen_func.denominator().factor()
        cyclic_orders = []
        for p, e in factors:
            found = False
            for a in [1, 2, 3, 4, 5, 6]:
                if p == cyclotomic_polynomial(a):
                    cyclic_orders.append(a)
                    found = True
            assert found
        return lcm(cyclic_orders)

    def l1_norm(gen_func):
        p = gen_func.numerator()
        q = gen_func.denominator()
        return int(p.norm(1)) + int(q.norm(1))

    df['gen_func_exp_dim'] = df.gen_func_denom.apply(x_minus_one_factor) - 1
    df['gen_func_deg'] = df.gen_func.apply(gen_func_deg)
    df['gen_func_period'] = df.gen_func.apply(period)
    df['gen_func_l1_norm'] = df.gen_func.apply(l1_norm)

    gcd_factor = df.gen_func.apply(lambda x:gcd(x.numerator().coefficients()))
    df['gen_func_norm'] = df.gen_func/gcd_factor

    assert all(df.gen_func_exp_dim==df.LW_dim)

    df['LW_num_verts'] = df.LW_verts_comp.apply(lambda x:sum(eval(x)))
    return df




def examine_dim_1(row):
    assert row.LW_dim == 1
    vertices = eval(row.vertex_genera).keys()
    all_faces = eval(row.all_faces)
    edges = []
    for face in all_faces:
        if face['dim'] == 1:
           edges.append(face['verts'])


    G = nx.Graph()
    G.add_nodes_from(vertices)
    G.add_edges_from(edges)

    assert nx.is_forest(G)
    return max(dict(G.degree).values())



def is_simplicial(row):
    all_faces = eval(row['all_faces'])
    for face in all_faces:
        if len(face['verts']) != face['dim'] + 1:
            return False
    return True


class PolyhedralFace(object):
    def __init__(self, dim, vertices):
        self.dim = dim
        self.vertices = frozenset(vertices)
        self.codim_1_faces = []
        self._barycentric_simplices = None

    def is_face_of(self, other):
        return self.vertices.issubset(other.vertices)

    def _cpt_string(self):
        return '|'.join(sorted(self.vertices))

    def __repr__(self):
        return 'PF(%d; %s)' % (self.dim, self._cpt_string())

    def barycentric_simplices(self):
        if self._barycentric_simplices is None:
            if self.dim == 0:
                ans = [list(self.vertices)]
            else:
                v = self._cpt_string()
                ans = []
                for face in self.codim_1_faces:
                    for simplex in face.barycentric_simplices():
                        ans.append(simplex + [v])
            self._barycentric_simplices = ans
        return self._barycentric_simplices


def triangulate_LW(row):
    raw_faces = eval(row['all_faces'])
    if is_simplicial(row):
        return SimplicialComplex([face['verts'] for face in raw_faces])
    dim = int(row['LW_dim'])
    faces = {d:list() for d in range(dim + 1)}
    for raw_face in raw_faces:
        d = raw_face['dim']
        F = PolyhedralFace(d, raw_face['verts'])
        faces[d].append(F)

    for d in range(dim):
        for A in faces[d]:
            for B in faces[d + 1]:
                if A.is_face_of(B):
                    B.codim_1_faces.append(A)

    all_simplices = []
    all_face_names = []
    for d in range(dim + 1):
        for A in faces[d]:
            all_simplices += A.barycentric_simplices()
            all_face_names.append(A._cpt_string())

    X = SimplicialComplex(all_simplices)
    assert set(X.vertices()) == set(all_face_names)
    return X



def test_simplicial_complex(row):
    X = triangulate_LW(row)
    one_vert_per_comp = [C[0] for C in X.graph().connected_components()]
    for v in one_vert_per_comp:
        Y = X.connected_component([v])
        print(Y)
        print('Checking homology...')
        assert Y.is_acyclic()
        print('Checking fundamental group...')
        if Y.dimension() > 1:
            assert Y.fundamental_group().order() == 1

        if Y.dimension() > 0:
            P = Y.face_poset()
            n = Y.dimension()
            for face in Y.faces()[n-1]:
                face = tuple(face)
                assert len(P.upper_covers(face)) <= 2

        print('Skipping shellability test')
        # assert Y.is_shellable()

#def triangulate_LW(row):
#    vertices = eval(row.vertex_genera).keys()
#    all_faces = eval(row.all_faces)
#    maximal_faces = eval(row.maximal_faces)
#    if is_simplicial(row):
#        return SimplicialComplex([face['verts'] for face in maximal_faces])


def check_table_4(df):
    names = ['K10n10', 'K14n11913', 't12766', 'K15n93515', 'K12n605',
             'K11n34', 'K14n1808', 'K12n214', 'K15n15582']
    count = 0
    for name in names:
        row = df[df.name==name].iloc[0]
        B = row.gen_func()
        print('%d  %s / %s   %d %d %s %d\n\n' % (
            row.LW_dim, B.numerator(), row.gen_func_denom,
            row.gen_func_period, row.gen_func_l1_norm,
            row['name'], sum(df.gen_func==B)))
        count += sum(df.gen_func==B)

    return count



def check_table_5(df):
    names = ['K15n138922', 'K15n27228', 'K15n86383', 'K15n139871',
             'K13n1795', 'K13n2458']

    count = 0
    for name in names:
        row = df[df.name==name].iloc[0]
        B = row.gen_func()
        print('%d  %s / %s   %d %d %s %d\n\n' % (
            row.LW_dim, B.numerator(), row.gen_func_denom,
            row.gen_func_period, row.gen_func_l1_norm,
            row['name'], sum(df.gen_func==B)))
        count += sum(df.gen_func==B)

    return count


def format_by_genus_count_table(df):
    to_show = df.by_genus.value_counts().head(10)
    c = 0
    for a_g in to_show.index:
        M = df.name[df.by_genus==a_g].iloc[0]
        c += to_show[a_g]
        print(a_g[1:-1] + ' &   & ' + M + ' & %d \\\\' % c)
    return c

ag_common = ['t09753', 't12198', 'K14n11913', 'K12n605', 'K11n73', 'K14n13645',
             'K11n42', 'K11n34', 'o9_37085', 'K15n93515']

ag_messy = """\
[8,14,46,89,224,305,674,905,1536,1955,3326,3771,6150,7019,9850,11611,16714,17767,25490,27415]
[8,16,54,98,264,318,806,984,1794,2098,3994,4074,7368,7632,11552,12976,20114,19396,30670,30550]
[8,11,36,69,176,221,552,662,1200,1448,2760,2504,5210,5065,7576,8402,14114,11989,21676,18758]
[12,21,61,109,261,320,721,880,1480,1762,3094,3115,5429,5666,8019,9086,13596,13059,20062,19841]
[8,9,28,53,132,179,410,529,958,1153,2096,2233,3976,4101,6424,6899,11012,10393,16938,16443]
[10,25,71,140,352,473,1058,1386,2389,2939,5152,5585,9422,10311,14887,17057,25304,25573,38238,39603]
[12,16,51,99,235,345,711,999,1649,2209,3551,4319,6593,7919,10971,13231,18275,20555,28063,31485]
[8,18,57,110,270,356,785,1013,1737,2092,3667,3942,6614,7134,10397,11710,17426,17422,26131,26891]
[12,34,110,216,532,708,1558,2018,3462,4176,7314,7876,13204,14256,20778,23404,34820,34832,52226,53766]
[12,30,109,231,549,861,1737,2511,4059,5643,8859,10941,16623,20229,27303,33729,46215,52455,71079,80271]
[10,21,73,143,385,513,1224,1605,2870,3542,6409,7010,12051,13231,19463,22436,33614,34307,51700,53862]
""".split()

ag_messy_names = ['K12n214', 'K12n210', 'K13n866', 'K13n3763', 'K13n1019', 'K15n15582',
                  'K15n15220', 'K15n23198', 'K13n3838', 'K15n33595', 'K13n2458']

def format_by_genus_messy(dg, chars=80):
    names = []
    for a_g in ag_messy:
        match = dg[dg.by_genus.apply(lambda x:x.startswith(a_g[:-1]))].iloc[0]        
        M = match['name']
        names.append(M)
        print(a_g[1:-1].replace(',', ', ') + ' & ' + M + ' \\\\ \\aglinesp')
    return names

def moebius_invert(coeffs):
    """
    The Dirichlet convolution of the input with the Moebius function
    """
    n = len(coeffs)
    return [sum([coeffs[d - 1] * moebius(m/d) for d in divisors(m)])
            for m in range(1, n + 1)]

def moebius_transform(coeffs):
    """
    The Dirichlet convolution of the input with 1
    """
    n = len(coeffs)
    return [sum([coeffs[d - 1] for d in divisors(m)])
            for m in range(1, n + 1)]

def gen_func_of_seq(coeffs):
    """
    Generating function for given coeffs [c1, c2, ... ] with no
    constant term, e.g. c1 x + c2 x^2 + ...
    """
    n_max = len(coeffs) + 1
    R = PowerSeriesRing(QQ, 'x', default_prec=n_max)
    x = R.gen()
    ans = sum([a * x**(n+1) for n, a in enumerate(coeffs)])
    return ans + R(0, n_max)

def lambert_transform(coeffs):
    n_max = len(coeffs) + 1
    R = PowerSeriesRing(QQ, 'x', default_prec=n_max)
    x = R.gen()
    ans = sum([a * x**(n+1) /(1 - x**(n + 1)) for n, a in enumerate(coeffs)])
    return ans + R(0, n_max)

def rational_function_size(f):
    a, b = f.numerator(), f.denominator()
    A, B = a.coefficients(), b.coefficients()
    return sum(c.global_height() for c in A + B)

def check_for_regularity(coeffs):
    F = gen_func_of_seq(coeffs)
    LA = lambert_transform(coeffs)
    R = F.parent()
    n = (len(coeffs) - 1)//2
    guess = LA.pade(n, n)
    # Giant jump in data here from 14 to 800
    if rational_function_size(guess) < 20:
        assert (LA - R(guess)).is_zero()
        return guess
    else:
        assert not (LA - R(guess)).is_zero()

def format_gen_func(f):
    a, b = f.numerator(), f.denominator()
    b = b.factor()
    ans = '(' + repr(a) + ')/'
    if len(b) == 1:
        ans += repr(b)
    else:
        ans += '(' + repr(b) + ')'
    assert f.parent()(ans) == f
    return ans

def add_regularity():
    dg = pd.read_csv('data_LW/extended_by_genus.csv')
    dg['likely_regular'] = 0
    dg['LA_gen_func'] = None

    for i, row in dg.iterrows():
        coeffs = eval(row.by_genus)
        guess = check_for_regularity(coeffs)
        if not guess is None:
            dg.loc[i, 'likely_regular'] = 1
            dg.loc[i, 'LA_gen_func'] = format_gen_func(guess)
    return dg
            
    
    
def func(row):
    coeffs = eval(row.by_genus)
    return check_for_regularity(coeffs)
    

def strange(A, B, g_max):
    assert A.genus==2 and B.genus==3
    for g in range(2, g_max+1):
        for a in range(g):
            k = g - a - 1
            if k % 2 == 0:
                b = k//2
                if a > 0 and b > 0:
                    F = a*A + b*B
                    if F.connected:
                        assert F.genus == g
                        print(g, a, b)
                    

table_13_exs = ['t09753', 't12198', 'K14n11913', 'K12n605', 'K11n73', 'K15n67261', 'K15n129923', 'K15n138922']

def table_13_info():
    R = PolynomialRing(QQ, 'x')
    F = R.fraction_field()
    
    def gen_func_deg(gen_func):
        n = gen_func.numerator().degree()
        d = gen_func.denominator().degree()
        assert n == d
        return n

    def period(gen_func):
        factors = gen_func.denominator().factor()
        cyclic_orders = []
        for p, e in factors:
            found = False
            for a in [1, 2, 3, 4, 5, 6]:
                if p == cyclotomic_polynomial(a):
                    cyclic_orders.append(a)
                    found = True
            assert found
        return lcm(cyclic_orders)

    def l1_norm(gen_func):
        p = gen_func.numerator()
        q = gen_func.denominator()
        return int(p.norm(1)) + int(q.norm(1))
    
    dg = pd.read_csv('data_LW/extended_by_genus.csv')
    for M in table_13_exs:
        row = dg[dg.name==M].iloc[0]
        gen = row.LA_gen_func
        f = F(gen)
        print(gen.replace('*', '') + ' & ' + repr(period(f)) + ' & ' + repr(l1_norm(f)) + ' & ' + M + ' \\\\')
        

        
        
