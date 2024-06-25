# is_b_had.py
import awkward as ak
import numpy as np
import numba as nb
from memory_profiler import profile
@nb.njit
def is_b_had(pdgId, builder):
    ''' numba implementation for a function that checks if a particle is a B hadron according to its pdgId
    Input: pdgId: a 2d awkward array of integers
    Output: a 2d awkward array of [True, False] if the particle is a B hadron or not
    '''
    for row in pdgId:
        builder.begin_list()
        for xii in row:
            x=str(np.abs(xii))
#             if len(x)>2 and (x[0]=='5' or (len(x)>4 and x[-3]=='5') ):
#                 print(xii)
            if len(x)>2 and (x[0]=='5' or (len(x)>4 and x[-3]=='5') ):
#                 res = 1
#                 print(xii)
                builder.integer(1)
            else:
#                 res = 0
                builder.integer(0)
                
#             if res == 1:
#                 builder.integer(xii)
#             else:
#                 builder.integer(0)
                
#             builder.end_list()
        builder.end_list()
    return builder

@nb.njit
def is_b_had_3d(pdgId, builder):
    ''' Same as is_b_had but for a 3d array
    '''
    for row in pdgId:
        builder.begin_list()
        for row2d in row:
#             print(row)
            builder.begin_list()
            if len(row2d)>0:
                for xii in row2d:
#                     print(xii)
                    x=str(np.abs(xii))
            #             if len(x)>2 and (x[0]=='5' or (len(x)>4 and x[-3]=='5') ):
            #                 print(xii)
                    if len(x)>2 and (x[0]=='5' or (len(x)>4 and x[-3]=='5') ):
            #                 res = 1
            #                 print(xii)
                        builder.integer(1)
                    else:
                        builder.integer(0)
            builder.end_list()
        builder.end_list()
    return builder

def vectorial_is_b_bar(x):
    ''' A much simpler implmentation of the function above
    '''
    x = np.abs(x)
    x_mod_100 = x//100
    x_mod_1000 = x//1000
    res = np.zeros_like(x)
    test1 = x_mod_100 == 5
    test2 = x_mod_1000 == 5
    test3 = np.logical_and((x_mod_100>100), (x_mod_100 % 10)==5)
    return test1 | test2 | test3

def vectorial_is_c_car(x):
    ''' A much simpler implmentation of the function above
    '''
    x = np.abs(x)
    x_mod_100 = x//100
    x_mod_1000 = x//1000
    res = np.zeros_like(x)
    test1 = x_mod_100 == 4
    test2 = x_mod_1000 == 4
    test3 = np.logical_and((x_mod_100>100), (x_mod_100 % 10)==4)
    return test1 | test2 | test3

# @profile
def find_gluon_split_jets(genpart, b_jets, dr_cut=0.4, quark_type='b'):
    ''' A function that finds the gluon splitting jets by checking if there are at least two b/c hadrons within a radius dr_cut from the jet.
    quark_type: 'b' or 'c'; type of the heavy quark searched for.
    '''
    if quark_type=='b':
        check_hadron = vectorial_is_b_bar
    elif quark_type=='c':
        check_hadron = vectorial_is_c_car
    else:
        raise ValueError(f"Wrong value for quark type. Given value: {quark_type}. Available: b or c.")

    b_had_bool = check_hadron(genpart.pdgId)
    bhads = genpart[b_had_bool[:]==1]
    # print("bhads", bhads[:3])
    daughter_is_b = check_hadron(bhads.children.pdgId)
    # print("daughter_is_b", daughter_is_b[:3])
    no_b_daughter = ak.sum(daughter_is_b,axis=2)==0
    final_bhads = bhads[no_b_daughter]

    drs = b_jets.metric_table(final_bhads)
    is_from_gluon_splitting = np.sum(drs<dr_cut,axis=2)>=2
    is_not_from_gluon_splitting = np.sum(drs<dr_cut,axis=2)==1
    return is_from_gluon_splitting, is_not_from_gluon_splitting
    # return np.ones_like(b_jets.pt), np.zeros_like(b_jets.pt)

# # @profile
# def find_gluon_split_jets(genpart, b_jets, dr_cut=0.4):
#     # # print("genpart", genpart[:3])
#     # # print("b_jets", b_jets[:3])
#     # aa = genpart.pdgId
#     b_had_bool = is_b_had(genpart.pdgId, ak.ArrayBuilder())
#     bhads = genpart[b_had_bool[:]==1]
#     print("bhads", bhads[:3])
#     daughter_is_b = is_b_had_3d(bhads.children.pdgId, ak.ArrayBuilder())
#     print("daughter_is_b", daughter_is_b[:3])
#     no_b_daughter = ak.sum(daughter_is_b,axis=2)==0
#     final_bhads = bhads[no_b_daughter]

#     drs = b_jets.metric_table(final_bhads)
#     is_from_gluon_splitting = np.sum(drs<dr_cut,axis=2)==2
#     is_not_from_gluon_splitting = np.sum(drs<dr_cut,axis=2)==1
#     return is_from_gluon_splitting, is_not_from_gluon_splitting
#     # return np.ones_like(b_jets.pt), np.zeros_like(b_jets.pt)
