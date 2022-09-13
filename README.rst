==========================================
Counting essential surfaces in 3-manifolds
==========================================

This code and data accompanies the paper of the same name by Nathan
Dunfield, Stavros Garoufalidis, and Hyam Rubinstein.

Data
====

The "data" folder contains the raw data summarized in Sections 7 and 8
of the paper. These are all compressed "csv" files which can be opened
in any spreadsheet program (e.g. Excel, Libre Office) after being
uncompressed. For more programmatic exploration, we recommend Python's
Pandas data-wrangling module. The individual files are as follows:

1. The most interesting are "very_large_combined.csv" and
   "extended_by_genus.csv". The former describes the 4,330 very large
   examples discussed in Section 7 and the latter the very detailed genus
   counts studied in Section 8.

2. "cusped_census_all.csv" and "knots_large.csv" include the barely
   large, and, in the first case, the small manifolds in the samples.

3. "cusped_nonnormal_iso.csv" and "knots_nonnormal_iso.csv" give
   examples of triangulations where there are nonnormal isotopies of
   lw-surfaces.

4. "surface_counts_stress.csv" includes multiple triangulations of
   each manifold and was used to help validate the correctness of the
   algorithm as described in Section 7.7.

Here are descriptions of some of the key columns in these files:

a. tri_used: the triangulation isosig of the triangulation T for which
   LW_T was computed.  Specifically, these were loaded into Regina via::

     regina.Triangulation3(snappy.Triangulation(isosig)._to_string())

   which means that the triangulation uses was oriented and not just
   orientable.

b. gen_func: The rational form P(x)/Q(x) of the generating function
   B_M(x) from Theorem 1.3.

c. by_genus: The sequence a_M(g) starting at g = 2 and going up to
   g = 21. The file "extended_by_genus.csv" goes further (typically)
   and includes whether a_M(g) appears regular and, if so, what our
   conjectured generating function LA_M(x) is.

d. LW_dim, LW_euler, LW_num_max_faces, etc: record properties of LW_T.

e. vertex_surfaces, max_faces, all_faces: This is a detailed
   description of the complex LW_T. We provide code for quickly
   recovering LW_T from this condensed form, see below.
  
f. LW_euler_bound: 2 - 2 g_0 where g_0 is as Algorithm 6.12.

g. blah_time: time, in seconds, needed to compute the contents of
   column "blah", in 2020 on not particularly new hardware (typically
   single-threaded).

h. small, regina_small, nonnormal_iso: Records the indicated property
   where 1 means it has the property, -1 means it does not, and 0
   means unknown.


Code quickstart
===============

The included "nscomplex" Python module requires SageMath with SnapPy,
Regina, and the optional Normaliz module installed. Tested with
SageMath 8.9 and 9.1, SnapPy 2.7 and 2.8, and Regina 5.1. For working
with the data files, Pandas is also recommended. An easy way to get
all these moving parts working together is to use the Kitchen Sink::

  https://snappy.computop.org/installing.html#kitchen-sink

To test if everything is working, from within this directory type::

  sage -python -m nscomplex.test

It should print various messages ending with (in < 1 minute)::

  All doctests:
     0 failures out of 165 tests.
  
Here is a sample session that assumes you have started Sage from
within this directory::

  sage: import nscomplex
  sage: CS = nscomplex.ConnectedSurfaces('t12071', euler_bound=-6)
  sage: len(CS.normal), len(CS.almost_normal), len(CS.tubed), len(CS.least_weight)
  (44, 28, 175, 30)
  sage: LW = CS.essential_faces_of_normal_polytope(); LW
  FacesComplex(dim=2, num max=6, num vert=10, euler=1)
  sage: LW.ehrhart_series()     # This is B_M(x)
  (-x^3 + x^2 - 10*x)/(x^3 - 3*x^2 + 3*x - 1)
  sage: LW.num_of_genus(2, 11)  # a_M(g) for g= 2,...,10
  [10, 8, 12, 12, 20, 12, 28, 20, 28]

If you want to be able to use the "nscomplex" module from other
directories, you can install it into Sage as follows::

  sage -pip install .


Using stored LW_T
=================

Here is an example of rehydrating a stored LW_T, which took 70 hours
to compute from scratch::

  sage: import nscomplex
  sage: import pandas as pd
  sage: df = pd.read_csv('data/very_large_combined.csv')
  sage: i = df.ess_face_time.idxmax()
  sage: row = df.loc[i]
  sage: row
  id                                                               12234
  name                                                          K14n2035
  volume                                                         15.6716
  tets                                                                17
  small                                                               -1
  regina_small                                                       NaN
  nonnormal_iso                                                        0
  LW_dim                                                               4
  LW_euler                                                             1
  LW_comps                                                             1
  LW_dim_comps                                                       [4]
  LW_euler_comps                                                     [1]
  LW_inhom_comp                                                        0
  LW_max_vert_face                                                    10
  LW_verts_comp                                                     [22]
  LW_num_max_faces                                                    12
  tri_used             rLLvALAzwAQQccdhfhjimnlpmqopoqqdctagaaaqkaotafowb
  ess_face_time                                                   256528
  regina_small_time                                                  NaN
  gen_func_time                                                  4.26226
  by_genus_time                                                   4423.6
  vertex_surfaces      {"N196":[0,1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0...
  vertex_genera        {"N196":4,"N49":3,"N3":2,"N44":3,"N29":2,"N22"...
  max_faces            [{"dim":4,"verts":["N13","N15","N22","N24","N8...
  all_faces            [{"dim":0,"verts":["N50"]},{"dim":0,"verts":["...
  LW_euler_bound                                                     -16
  gen_func             (-x^6 + 4*x^5 - 6*x^4 - 3*x^2 - 12*x)/(x^6 - 4...
  by_genus             [12,21,61,109,261,320,721,880,1480,1762,3094,3...
  by_genus_naive                          [12,21,61,109,261,320,721,880]
  Name: 3441, dtype: object

  sage: surfaces, LW = nscomplex.reconstruct_faces(row)
  sage: LW
  FacesComplex(dim=4, num max=12, num vert=22, euler=1)
  sage: LW.maximal
  [AFace(d=4, [N13, N15, N22, N24, N87, N125]),
   AFace(d=4, [N2, N3, N17, N49, N50]),
   AFace(d=4, [N2, N13, N24, N44, N87]),
   AFace(d=4, [N2, N3, N17, N32, N49, N196]),
   AFace(d=4, [N2, N3, N6, N49, N196]),
   AFace(d=4, [N13, N22, N24, N44, N87]),
   AFace(d=4, [N2, N6, N13, N15, N16, N24, N37, N87, N125, N126]),
   AFace(d=4, [N15, N16, N24, N26, N37, N126]),
   AFace(d=4, [N2, N3, N17, N32, N56]),
   AFace(d=4, [N2, N3, N29, N32, N56]),
   AFace(d=4, [N2, N3, N6, N11, N49, N50]),
   AFace(d=4, [N15, N22, N24, N26, N125, N126])]
  sage: LW.ehrhart_series()
  (-x^6 + 4*x^5 - 6*x^4 - 3*x^2 - 12*x)/(x^6 - 4*x^5 + 5*x^4 - 5*x^2 + 4*x - 1)


Code details
============

In the "nscomplex" directory, the main files are as follows: 

* regina_util.py: As the names suggests, helper code for getting
  information in and out of Regina via the latter's Python interface.

* enumerate_surfaces.py: Uses Regina to enumerate all normal surfaces
  and almost normal surfaces with octagons with Euler characteristic
  bounded below, starting with the lists of fundamental surfaces.

* surfaces.py: A key file, this provides convenience wrappers
  NormalSurface and AlmostNormalSurface for Reginas normal surfaces
  and almost normal surfaces with octagons. Also provides TubedSurface
  for almost normal surfaces with tubes (not present in Regina). The
  function "connected_surfaces_to_euler_char" does as the name suggests,
  returning the answer as NormalSurfaces, AlmostNormalSurfaces, and
  TubedSurfaces. The enumeration of the first two is handled by
  "enumerate_surfaces.py" with the equivalence classes of tubes worked
  out by "skelata.py".
  
* skeleta.py: Implements the tightening of almost normal surfaces as
  well as equivalence classes of tubed surfaces associated to a given
  normal surface.

* faces.py: Working with faces of the normal surface polytope P_T and
  subcomplexes of it such as LW_T.

* connected_surface.py: The top level, this file defines the
  ConnectedSurface class which implements the algorithms in Section 6
  of the paper.

In the "scripts" directory is some of the code used to generate and
analyze the data, included as examples of the nscomplex code in
action.  Not everything in these files will work for you as there are
some private dependencies.

License
=======

All code and data herein is hereby released into the public domain by
its above-named authors, as per CC0::

  https://creativecommons.org/publicdomain/zero/1.0/legalcode


