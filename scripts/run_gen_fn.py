#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j
#

import nscomplex
import regina
import snappy
import json, time
# import taskdb2.worker
from sage.all import vector
import generate_tris
from nscomplex import reconstruct_faces

def find_essential_faces(manifold):
    CS = nscomplex.EssentialSurfaces(manifold)
    try:
        faces = CS.essential_faces_of_normal_polytope()
    except nscomplex.MoreSurfacesNeedToFindFaces:
        euler_bound = CS.euler_bound
        print('Failed at euler_bound: %d ' % CS.euler_bound)
        while True:
            euler_bound += -2
            CS = nscomplex.EssentialSurfaces(manifold, euler_bound=euler_bound)
            try:
                faces = CS.essential_faces_of_normal_polytope()
                break
            except nscomplex.MoreSurfacesNeedToFindFaces:
                print('Failed at euler_bound: %d ' % CS.euler_bound)

    return CS, faces

def summarize_essential_faces(faces):
    vertex_surfaces = set.union(*[set(face.vertex_surfaces)
                                  for face in faces.maximal])
    vertex_surfaces = sorted(vertex_surfaces, key=lambda S:S.index)

    def format_face(face):
        return {'dim':face.dim, 'verts':[repr(S) for S in face.vertex_surfaces]}

    data = {'LW_dim':max(face.dim for face in faces.maximal),
            'num_max_faces':len(faces.maximal),
            'vertex_surfaces':{repr(S):S.quad_vector for S in vertex_surfaces},
            'vertex_genera':{repr(S):S.genus for S in vertex_surfaces},
            'max_faces':[format_face(face) for face in faces.maximal],
            'all_faces':[format_face(face) for face in faces.faces],
            'LW_euler':faces.euler}

    return data

def add_essential_face_info(task):
    start = time.time()
    M0 = snappy.Triangulation(task['name'])
    isosig = M0.triangulation_isosig(decorated=False)
    CS, faces = find_essential_faces(isosig)
    if len(faces) == 0:
        task['small'] = 1
    else:
        task['small'] = -1
        task['LW_euler_bound'] = CS.euler_bound
        summary = summarize_essential_faces(faces)
        for col in ['LW_dim', 'num_max_faces', 'LW_euler']:
            task[col] = summary[col]
        for col in ['vertex_genera','vertex_surfaces']:
            task[col] = json.dumps(summary[col]).replace(' ', '')
        for col in ['max_faces', 'all_faces']:
            task[col] = json.dumps(summary[col]).replace(' ', '')
        task['nonnormal_iso'] = 1 if CS.has_isotopy_of_least_weight() else 0
        naive_counts = CS.num_incompressible_by_genus(concise=True)
        task['by_genus_naive'] = repr(naive_counts).replace(' ', '')
    run_time = time.time() - start

    task['ess_face_time'] = run_time
    task['tri_used'] = isosig
    task['done'] = True

def add_essential_face_info_try_harder(task):
    start = time.time()
    M0 = snappy.Triangulation(task['name'])
    tris = generate_tris.many_triangulations(M0, 60)
    isosigs = [T.triangulation_isosig(decorated=False) for T in tris]
    isosigs.remove(M0.triangulation_isosig(decorated=False))
    for isosig in isosigs[:15]:
        print('   Starting ' + isosig + ' for ' + task['name'])
        CS, faces = find_essential_faces(isosig)
        assert len(faces) > 0
        if not CS.has_isotopy_of_least_weight():
            if len(faces) == 0:
                task['small'] = 1
            else:
                task['small'] = -1
                task['LW_euler_bound'] = CS.euler_bound
                summary = summarize_essential_faces(faces)
                for col in ['LW_dim', 'num_max_faces', 'LW_euler']:
                    task[col] = summary[col]
                for col in ['vertex_genera','vertex_surfaces']:
                    task[col] = json.dumps(summary[col]).replace(' ', '')
                for col in ['max_faces', 'all_faces']:
                    task[col] = json.dumps(summary[col]).replace(' ', '')
                task['nonnormal_iso'] = 1 if CS.has_isotopy_of_least_weight() else 0
                naive_counts = CS.num_incompressible_by_genus(concise=True)
                task['by_genus_naive'] = repr(naive_counts).replace(' ', '')
            run_time = time.time() - start

            task['ess_face_time'] = run_time
            task['tri_used'] = isosig
            task['tets'] = snappy.Triangulation(isosig).num_tetrahedra()
            task['done'] = True
            return

def add_gen_func(task):
    """
    Needs columns:

    ['nonnormal_iso', 'vertex_genera', 'vertex_surfaces', 'all_faces', 'tri_used']
    """
    assert task['nonnormal_iso'] == 0
    vertex_surfaces, all_faces = reconstruct_faces(task)
    start = time.time()
    task['gen_func'] = repr(sum(face.ehrhart_series_of_interior() for face in all_faces))
    run_time = time.time() - start
    task['gen_func_time'] = run_time
    task['done'] = True

def add_count_by_genus(task):
    """
    Needs columns:

    ['nonnormal_iso', 'vertex_genera', 'vertex_surfaces', 'all_faces', 'tri_used']
    """
    assert task['nonnormal_iso'] == 0
    vertex_surfaces, all_faces = reconstruct_faces(task)
    start = time.time()

    def num_of_genus(g):
        return sum([face.num_of_genus_in_interior(g) for face in all_faces])

    task['by_genus'] = repr([num_of_genus(g) for g in range(2, 22)]).replace(' ', '')
    run_time = time.time() - start
    task['by_genus_time'] = run_time
    task['done'] = True

def add_count_by_genus_50(task):
    """
    Needs columns:

    ['nonnormal_iso', 'vertex_genera', 'vertex_surfaces', 'all_faces', 'tri_used']
    """
    assert task['nonnormal_iso'] == 0
    vertex_surfaces, all_faces = reconstruct_faces(task)
    start = time.time()

    def num_of_genus(g):
        return sum([face.num_of_genus_in_interior(g) for face in all_faces])

    task['by_genus_50'] = repr([num_of_genus(g) for g in range(2, 52)]).replace(' ', '')
    run_time = time.time() - start
    task['by_genus_50_time'] = run_time
    task['done'] = True


def add_count_by_genus_100(task):
    """
    Needs columns:

    ['nonnormal_iso', 'vertex_genera', 'vertex_surfaces', 'all_faces', 'tri_used']
    """
    assert task['nonnormal_iso'] == 0
    vertex_surfaces, all_faces = reconstruct_faces(task)
    start = time.time()

    def num_of_genus(g):
        return sum([face.num_of_genus_in_interior(g) for face in all_faces])

    task['by_genus_100'] = repr([num_of_genus(g) for g in range(2, 102)]).replace(' ', '')
    run_time = time.time() - start
    task['by_genus_100_time'] = run_time
    task['done'] = True


def add_count_by_genus_200(task):
    """
    Needs columns:

    ['nonnormal_iso', 'vertex_genera', 'vertex_surfaces', 'all_faces', 'tri_used']
    """
    assert task['nonnormal_iso'] == 0
    vertex_surfaces, all_faces = reconstruct_faces(task)
    start = time.time()

    def num_of_genus(g):
        return sum([face.num_of_genus_in_interior(g) for face in all_faces])

    task['by_genus_200'] = repr([num_of_genus(g) for g in range(101, 202)]).replace(' ', '')
    run_time = time.time() - start
    task['by_genus_200_time'] = run_time
    task['done'] = True
    

def add_regina_small(task):
    start = time.time()
    T = nscomplex.regina_util.as_regina(str(task['name']))
    start = time.time()
    small = nscomplex.regina_util.is_small(T)
    task['regina_small'] = 1 if small else -1
    run_time = time.time() - start
    task['regina_small_time'] = run_time
    task['done'] = True

def add_LW_dull(task):
    import analysis
    analysis.test_simplicial_complex(task)
    task['LW_dull'] = 1
    task['done'] = True
    #except:
    #    print('FAILED')
    #    assert False
    
    
        
    
# class AdmissibleFaceRestored(AdmissibleFace):
#     def __init__(self, triangulation, face_data):
#         T = triangulation
#         vert_surfaces = []
#         sum_vector = 0
#         for label, quad_vec in face_data['verts'].items():
#             assert label[0] == 'N'
#             sum_vector += vector(quad_vec)
#             R = regina.NormalSurface(T, regina.NS_QUAD_CLOSED, list(quad_vec))
#             S = surfaces.Surface(R, int(label[:1]))
#             vert_surfaces.append(R)
#         AdmisibleFace.__init__(self, face_data['dim'], zero_set(sum_vector), vert_surfaces)

#     pass



#task = {'name':'K10a78'}
task1 = {'name':'K10n10'}
task2 = {'name':'kLLLMLQkcdfefhhijjjifeewvjqxok'}
task3 = {'name':'m004'}

#taskdb2.worker.run_function('surface_counts_cusped_new', 'task_ess_faces', add_essential_face_info)
#taskdb2.worker.run_function('surface_counts_knots_new', 'task_gen_func', add_gen_func)
#taskdb2.worker.run_function('surface_counts_knots_new', 'task_by_genus', add_count_by_genus)
#taskdb2.worker.run_function('surface_counts_cusped', 'task_regina', add_regina_small)
#taskdb2.worker.run_function('surface_counts_knots_new', 'task_try_harder', add_essential_face_info_try_harder)
#taskdb2.worker.run_function('surface_counts_deep', 'task_100', add_count_by_genus_100)
#taskdb2.worker.run_function('surface_counts_deep', 'task_200', add_count_by_genus_200)
#taskdb2.worker.run_function('surface_counts_knots_new', 'task_dull', add_LW_dull)
