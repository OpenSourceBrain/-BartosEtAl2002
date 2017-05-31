execfile('makeVidaCell.py')
makeCell('VidaCell','Vida.cell_morphology.nml')
execfile('makeVidaNetwork.py')
# to run simulation in commandline:
# jnml LEMS_Sim_VidaTest.xml -neuron -run -nogui

# should produce a file called v_BC.dat:
# format: time(s), v_1(V), v_2(V), ..., v_N(V)
# v_i=soma potential of i-th cell in V
