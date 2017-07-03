'''
Code to generate networks of Wang-Buszaki-type interneurons as described in 

"Fast synaptic inhibition promotes synchronized gamma oscillations in hippocampal interneuron networks"
Marlene Bartos, Imre Vida, Michael Frotscher, Axel Meyer, Hannah Monyer, Joerg R.P. Geiger, and Peter Jonas
PNAS 99(20):13222-13227 (2002) DOI 10.1073/pnas.192233099

"Shunting inhibition improves robustness of gamma oscillations in hippocampal interneuron networks by homogenizing firing rates"
Imre Vida, Marlene Bartos, and Peter Jonas
Neuron 49, 107-117 (2006) DOI 10.1015/j.neuron.2005.11.036

Adpated to existing NEURON model published in ModelDB:
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=21329

Python code adapted and modified from existing code in 
https://github.com/OpenSourceBrain/OpenCortex

Author: Birgit Kriener, April 2017
'''


from neuroml import __version__
from neuroml import NeuroMLDocument
from neuroml import Network
from neuroml import Population
from neuroml import Location
from neuroml import Instance
from neuroml import Projection
from neuroml import Property
from neuroml import PulseGenerator
from neuroml import GapJunction
from neuroml import Connection
from neuroml import ConnectionWD
from neuroml import IncludeType
from neuroml import InputList
from neuroml import Input
from neuroml import PoissonFiringSynapse
from neuroml import ExpTwoSynapse
from neuroml import ElectricalConnectionInstance
from neuroml import ElectricalProjection
from neuroml import VoltageClamp
import neuroml.writers as writers
from pyneuroml import pynml
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
from random import random
from random import uniform
from random import gauss
from random import seed
from math import exp



def add_gap_junction_synapse(nml_doc, id, conductance):
    
    """
    Adds a <gapJunction> element to the document. See the definition of the 
    behaviour of this here: https://www.neuroml.org/NeuroML2CoreTypes/Synapses.html#gapJunction
    
    Returns the class created.
    """
    
    syn0 = GapJunction(id=id, conductance=conductance)

    nml_doc.gap_junctions.append(syn0)

    return syn0


def add_elect_connection(projection, id, presynaptic_population, pre_component, pre_cell_id, pre_seg_id, postsynaptic_population, post_component, post_cell_id, post_seg_id, gap_junction_id, pre_fraction=0.5, post_fraction=0.5):

    """
    Add a single electrical connection (via a gap junction) to a projection between `presynaptic_population` and `postsynaptic_population`
    """

    connection = ElectricalConnectionInstance(id=id, \
                                              pre_cell="../%s/%i/%s" % (presynaptic_population, pre_cell_id, pre_component), \
                                              post_cell="../%s/%i/%s" % (postsynaptic_population, post_cell_id, post_component), \
                                              synapse=gap_junction_id, \
                                              pre_segment=pre_seg_id, \
                                              post_segment=post_seg_id, \
                                              pre_fraction_along=pre_fraction, \
                                              post_fraction_along=post_fraction)

    projection.electrical_connection_instances.append(connection)


def add_probabilistic_ring_gapjunx(net, presynaptic_population, pre_component, postsynaptic_population, post_component, prefix, gap_junction_id, numCells, numNeighbors_GJ):

    if numCells==0:
        return None
        
    proj = ElectricalProjection(id="%s_%s_%s"%(prefix,presynaptic_population, postsynaptic_population), presynaptic_population=presynaptic_population, postsynaptic_population=postsynaptic_population)

    count = 0

    for i in range(0, numCells):
        for j in range(-numNeighbors_GJ,numNeighbors_GJ+1):
            nn=(i+j)%numCells
            if (not i==nn) and (random()<0.5):
                add_elect_connection(proj, count, presynaptic_population, pre_component,i, 0, postsynaptic_population, post_component, nn, 0, gap_junction_id, pre_fraction=0.5, post_fraction=0.5)
                count+=1
    net.electrical_projections.append(proj)                
    return proj


def add_connectionWD(projection, id, pre_pop, pre_component, pre_cell_id, pre_seg_id, post_pop, post_component, post_cell_id, post_seg_id, weight, delay):

    connection = ConnectionWD(id=id, \
                            pre_cell_id="../%s/%i/%s"%(pre_pop, pre_cell_id, pre_component), \
                            pre_segment_id=pre_seg_id, \
                            pre_fraction_along=0.5,
                            post_cell_id="../%s/%i/%s"%(post_pop, post_cell_id, post_component), \
                            post_segment_id=post_seg_id,
                            post_fraction_along=0.5, delay=delay, weight=weight)

    projection.connection_wds.append(connection)
    

def add_probabilistic_ring_projection_delay(net, presynaptic_population, pre_component, postsynaptic_population, post_component, prefix, synapse, numCells, numNeighbors, weight):
    ### Gauss profile used in Vida et al (2006): f(d) = 9.9736*exp(d**2/-1152), i.e. footprint of sigma=24
    ### then compared to uniform random number in [0,10]; here we use probability profile f(d) in [0,1]
    def connprob(x,sig):
        return  exp(-x**2/2./sig**2)

    sigma = 24*numCells/200.  # chosen sigma to get ~57 neighbors assuming numCells=200
    d0    = 0.5               # constant synaptic delay in ms
    
    if numCells==0:
        return None
        
    proj = Projection(id="%s_%s_%s"%(prefix,presynaptic_population, postsynaptic_population), presynaptic_population=presynaptic_population, postsynaptic_population=postsynaptic_population, synapse=synapse)

    count = 0

    for i in range(0, numCells):
        for j in range(-numNeighbors,numNeighbors+1):
            nn=(i+j)%numCells
            if (not i==nn) and (random()<connprob(abs(j),sigma)):
                dring = abs(j)*0.2 # distance-dependent conduction delay in ms;
                                   # assuming pairwise distance=50um between 200 neurons
                                   # on the ring and conduction velocity of 0.25 m/s
                add_connectionWD(proj, count, presynaptic_population, pre_component, i, 0, postsynaptic_population, post_component, nn, 0, weight, str(d0+dring)+" ms")
                count+=1
                    
    net.projections.append(proj)
    print "connprob = ",count/(2.*numNeighbors*numCells)
    return proj
    


def add_population_in_rectangular_region(net, pop_id, cell_id, size, x_min, y_min, z_min, x_size, y_size, z_size, color=None):
    
        pop = Population(id=pop_id, component=cell_id, type="populationList", size=size)
        if color is not None:
            pop.properties.append(Property("color",color))
        net.populations.append(pop)

        for i in range(0, size) :
                index = i
                inst = Instance(id=index)
                pop.instances.append(inst)
                inst.location = Location(x=str(x_min +(x_size)*random()), y=str(y_min +(y_size)*random()), z=str(z_min+(z_size)*random()))

                

def generate_BC_cell_layer(network_id,
                        x_size,        # um
                        y_size,        # um
                        z_size,        # um
                        bc_group_component,
                        bc_syn_weight, # in nS
                        bc_syn_erev,   # in mV
                        bc_gj_weight,  # in nS
                        I_mu,          # mean input drive in [nA]
                        I_cv,          # std of input drive in terms of I_mu across neurons
                        numCells_bc = 200,
                        numNeighbors_bc = 50,
                        activate = False,
                        gapjunx = True,
                        validate = True,
                        random_seed = 12345,
                        generate_lems_simulation = False,
                        duration = 500, # in [ms]
                        dt = 0.01):
    
    seed(random_seed)

    nml_doc = NeuroMLDocument(id=network_id)

    net = Network(id = network_id)
                  
    net.notes = "Network generated using libNeuroML v%s"%__version__
    nml_doc.networks.append(net)

    if numCells_bc>0:
        nml_doc.includes.append(IncludeType(href='%s.cell.nml'%bc_group_component))

    ### The names of the groups/populations 
    bc_group = "BasketCells"

    ### Generate basket cells 
    if numCells_bc>0:
        add_population_in_rectangular_region(net, bc_group, bc_group_component, numCells_bc, 0, 0, 0, x_size, y_size, z_size, color="0 0 1")
        
    ### Connect cells in randomized ring structure
    ### Exp2Syn was used in Neuron code on modelDB
    ### If network shall be activated only after 150 ms; note:
    ### weight is set in expTwoSynapseAct.nml and cannot be set in call of projection routine!!
    if activate:
        nml_doc.includes.append(IncludeType('expTwoSynapseAct.nml'))
        add_probabilistic_ring_projection_delay(net, bc_group, bc_group_component, bc_group, bc_group_component, 'NetConn', 'expTwoSynapseAct', numCells_bc, numNeighbors_bc, 1)
    ### If network activated right away:    
    else:
        bc_syn = ExpTwoSynapse(id="bc_syn", gbase="1nS", erev=str(bc_syn_erev)+"mV", tau_rise="0.16 ms", tau_decay="1.8 ms")
        nml_doc.exp_two_synapses.append(bc_syn)
        add_probabilistic_ring_projection_delay(net, bc_group, bc_group_component, bc_group, bc_group_component, 'NetConn', "bc_syn", numCells_bc, numNeighbors_bc, bc_syn_weight)
 
    ### Add gap junctions randomly to 8 closes neighbors with prob=1/2    
    if gapjunx:
        gj_syn = add_gap_junction_synapse(nml_doc, id="gj_syn", conductance=str(bc_gj_weight)+"nS")
        add_probabilistic_ring_gapjunx(net, bc_group, bc_group_component, bc_group, bc_group_component, 'GapJunx', "gj_syn", numCells_bc, 4)

        
    ### Randomize initial potentials via voltage clamp
    vc_dur = 2  # ms
    for i in range(0, numCells_bc):
        tmp = -70. + random()*20
        vc = VoltageClamp(id='VClamp%i'%i, delay='0ms', duration='%ims'%vc_dur, simple_series_resistance='1e6ohm', target_voltage='%imV'%tmp)
        
        nml_doc.voltage_clamps.append(vc)
        
        input_list = InputList(id='input_%i'%i, component='VClamp%i'%i, populations=bc_group)
        input = Input(id=i, target='../%s/%i/%s'%(bc_group, i, bc_group_component), destination='synapses')
        input_list.input.append(input)
        
        net.input_lists.append(input_list)
        
    ### Define stimulus; f_surf relates to area of soma compartment; here we use area=100um^2 
    ### as in the modelDB demo code
    f_surf = 1e-3
    ### Gauss-distribution of driving currents 
    I_amps = [gauss(I_mu, I_cv*I_mu)*f_surf for _ in range(0,numCells_bc)]
    ### Uniform distribution of onset of driving currents in [ms]
    stim_delays = [uniform(0.,50.) for _ in range(0,numCells_bc)] 
    cnt = 0
    for i in range(0,numCells_bc):
        stim = PulseGenerator(id="stim%i"%(i),
                              delay=str(stim_delays[i])+'ms',
                              duration= str(duration)+'ms',
                              amplitude=str(I_amps[i])+'nA')
        
        nml_doc.pulse_generators.append(stim)

        input_list = InputList(id="%s_input"%stim.id,
                               component=stim.id,
                               populations=bc_group)
            
        syn_input = Input(id=cnt,
                          target="../%s/%i/%s" % (bc_group, i, bc_group_component),
                          destination="synapses")
        input_list.input.append(syn_input)
        net.input_lists.append(input_list)
        cnt+=1
    

    #######   Write to file  ######    

    print("Saving to file...")
    nml_file = network_id+'.net.nml'
    writers.NeuroMLWriter.write(nml_doc, nml_file)

    print("Written network file to: "+nml_file)


    if validate:

        ###### Validate the NeuroML ######    

        from neuroml.utils import validate_neuroml2
        validate_neuroml2(nml_file) 
        
    if generate_lems_simulation:
        ### Create a LEMSSimulation to manage creation of LEMS file
        
        ls = LEMSSimulation("Sim_%s"%network_id, duration, dt)

        ### Point to network as target of simulation
        ls.assign_simulation_target(net.id)
        
        ### Include generated/existing NeuroML2 files
        ls.include_neuroml2_file('%s.cell.nml'%bc_group_component)
        ls.include_neuroml2_file(nml_file)
        

        ### Specify Displays and Output Files
        if numCells_bc>0:
            disp_bc = "display_bc"
            ls.create_display(disp_bc, "Voltages Basket Cells", "-80", "40")

            of_bc = 'Volts_file_bc'
            ls.create_output_file(of_bc, "v_"+fname+".dat")


            for i in range(numCells_bc):
                quantity = "%s/%i/%s/v"%(bc_group, i, bc_group_component)
                ls.add_line_to_display(disp_bc, "MF %i: Vm"%i, quantity, "1mV", pynml.get_next_hex_color())
                ls.add_column_to_output_file(of_bc, "v_%i"%i, quantity)

        ### Save to LEMS XML file
        lems_file_name = ls.save_to_file()
    else:
        
        ls = None
        
    print "-----------------------------------"
    
    return nml_doc, ls

    
if __name__ == "__main__":
    
    fname = "VidaTest"    # file name
    numCells_bc = 100     # network size         
    numNeighbors_bc = int(numCells_bc/4.) # size of maximal neighborhood
    bc_syn_weight = 0.1   # peak conductance of interneuron synapses in [nS];
                          # NOTE: in original gmax_spec in [mS/cm^2],
                          # f_surf = 1e-3 cm^2, so gmax = |gmax_spec| in [nS]
    bc_syn_erev = -75.    # synaptic reversal potential in mV;                   
    activate = True       # bool to set whether network only activated after delay of 150 ms
    if not activate:      # if activate True: in expTwoSynapseAct.nml manually
        bc_syn_weight = bc_syn_weight # set g_base=bc_syn_weight in [nS]
        bc_syn_erev = bc_syn_erev     # set erev=bc_syn_erev in [mV]
    gapjunx = True    
    bc_gj_weight = 0.01   # conductance of gap junctions in [nS]
    I_mu = 1.             # mean input drive in [nA]
    I_cv = 0.             # std of input drive in terms of I_mu across neurons    
    x_size,y_size,z_size = 3*[500] # size of region that cells are distributed in [mum]
    duration = 500        # duration of simulation in [ms]
    dt = 0.1              # time step [ms]
   
    generate_BC_cell_layer(fname,
                           x_size = x_size, 
                           y_size = y_size, 
                           z_size = z_size,
                           bc_group_component = "VidaCell",
                           numCells_bc = numCells_bc,
                           numNeighbors_bc = numNeighbors_bc,
                           bc_syn_weight = bc_syn_weight,
                           bc_syn_erev = bc_syn_erev,
                           bc_gj_weight = bc_syn_weight,
                           I_mu = I_mu,
                           I_cv = I_cv,
                           activate = activate,
                           gapjunx = gapjunx,
                           generate_lems_simulation = True,
                           validate = True,
                           duration = duration,
                           dt = dt)


