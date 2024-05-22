function main()
    sys_dict = generate_inputs("InputFiles/Au_Input.txt")
    sim_settings = define_simulation_settings(sys_dict)
    material_parameters = define_material_parameters(sys_dict)
    

end

main()