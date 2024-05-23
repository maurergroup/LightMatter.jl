function main()
    #generate_inputs needs some work for restructuring
    sys_dict = generate_inputs("InputFiles/Au_Input.txt") #Dictionary of parameters
    sim_settings = define_simulation_settings(sys_dict) #Struct with sim settings
    material_parameters = define_material_parameters(sys_dict) #Struct with constant material parameters
    sim_dimension = define_sim_dimensions(sys_dict)
    laser = define_laser_system()
    sys = build_system()
end

main()