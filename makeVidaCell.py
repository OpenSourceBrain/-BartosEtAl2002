execfile('parameters_Vida.py')
def makeCell(outcellid, morphfile):
    import parameters_Vida as pm
    import neuroml
    import neuroml.loaders as loaders
    import neuroml.writers as writers
    ### cell morphology
    fn = morphfile
    doc = loaders.NeuroMLLoader.load(fn)
    print("Loaded morphology file from: "+fn)

    cell = doc.cells[0]
    
    ### channels and properties

    channel_densities = []

    cd_pas_soma = neuroml.ChannelDensity(id="pas_chan_soma", segment_groups="all", ion="non_specific", ion_channel="pas", erev=str(pm.erev)+" mV", cond_density=str(pm.gl)+" mS_per_cm2")
    channel_densities.append(cd_pas_soma)
    ### Na: soma
    cd_na_soma = neuroml.ChannelDensity(id="na_chan_soma", segment_groups="all", ion="non_specific", ion_channel="Na_BC", erev=str(pm.ena_soma)+" mV", cond_density=str(pm.gna_soma)+" mS_per_cm2")
    channel_densities.append(cd_na_soma)
    ### K: soma
    cd_k_soma = neuroml.ChannelDensity(id="k_chan_soma", segment_groups="all", ion="non_specific", ion_channel="K_BC", erev=str(pm.ek_soma)+" mV", cond_density=str(pm.gk_soma)+" mS_per_cm2")
    channel_densities.append(cd_k_soma)
    
    ### membrane properties

    specific_capacitances = []

    specific_capacitances.append(neuroml.SpecificCapacitance(value=str(pm.Cm)+' uF_per_cm2', segment_groups='all'))
    init_memb_potentials = [neuroml.InitMembPotential(value=str(pm.Vinit)+" mV", segment_groups='all')]
    spike_threshold = [neuroml.SpikeThresh(value=str(pm.spike_thresh)+" mV", segment_groups='all')]
    
    membrane_properties = neuroml.MembraneProperties(
        channel_densities=channel_densities,
        specific_capacitances=specific_capacitances,
        init_memb_potentials=init_memb_potentials,
        spike_threshes=spike_threshold)

    ### intracellular properties

    resistivities = []
    resistivities.append(neuroml.Resistivity(
        value=str(pm.Ri)+" ohm_cm", segment_groups='all'))

    intracellular_properties = neuroml.IntracellularProperties(resistivities=resistivities)

    bp = neuroml.BiophysicalProperties(id="biophys",
                                       intracellular_properties=intracellular_properties,
                                       membrane_properties=membrane_properties)
                                   
    cell.biophysical_properties = bp

    cell.id = outcellid

    nml_doc2 = neuroml.NeuroMLDocument(id=cell.id)

    nml_doc2.includes.append(neuroml.IncludeType('pas.channel.nml')) 
    nml_doc2.includes.append(neuroml.IncludeType('Na_BC.channel.nml')) 
    nml_doc2.includes.append(neuroml.IncludeType('K_BC.channel.nml'))

    nml_doc2.cells.append(cell)

    nml_file = cell.id+'.cell.nml'

    writers.NeuroMLWriter.write(nml_doc2,nml_file)
    
    print("Saved modified morphology file to: "+nml_file)
                                   
                                   
    ###### Validate the NeuroML ######    

    from neuroml.utils import validate_neuroml2
    
    validate_neuroml2(nml_file)



