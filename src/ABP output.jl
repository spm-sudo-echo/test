# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
# VARIABLES: Destination folder path and filename

using DrWatson
@quickactivate "active-brownian-particles"
using Plots,Distances,NaNStatistics,CSV, DataFrames
using FileIO, JLD2
using Dates

include(projectdir("src","ABP main.jl"))
include(projectdir("src","ABP file.jl"))
include(projectdir("src", projectdir("src", "measure", "ABP analysis.jl")))
include(projectdir("src","ABP multifolder.jl"))
# include(projectdir("src", "measure","ABP SD.jl"))


# -----------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL MAIN FUNCTION
# We plot the set of particles considering the correction of hard spheres

L = 100.0 	# μm box length
R = 2.0		# μm particle radius
v = 1.0 	# μm/s particle velocity
a=L/2
b=L/4
pf_factor = (R^2)
#pf_factor = (R^2)/(a*b)
DT, DR = diffusion_coeff(R).*[1e12, 1]
packing_fraction = 0.1
Np = round(Int,packing_fraction*L^2/(2R^2))  #Np is the number of particles in my set and I choose it random?
Nt = 1000   # Nt is the number of steps 


#-------------------------------------------------------------------------------------------------------------------
# RUN LOOP : multiple parameters sets and multiple runs for each (the output is saved in data/sims)
# parameters in the file must be a dict with list of values for each parameter 

function run(;param_file_name::String, nb_runs::Integer=1, wall_condition::String="periodic")
    # extract parameters instances from given file
    parameter_instances = load(datadir("parameters", param_file_name))
    # test if all parameter instances are complete (i.e. same number of parameters in each), return otherwise
    allequal(length.(values(parameter_instances))) || return

    # for each parameter instance 
    for parameters in parameter_instances
        # time signature to differentiate
        datestamp = Dates.format(now(),"YYYY-mm-dd_HH:MM:SS")  
        simulation_folder_path = datadir("sims", savename(datestamp, parameters; digits=3))
        mkdir(simulation_folder_path)
        # for each run
        # for i_run in 1:nb_runs
        #     file_path = joinpath(simulation_folder_path, "run_$i_run.jld2")
        #     mouvement_data = multiparticleE(;parameters...)
        #     @save file_path mouvement_data
        # end
    end
end

Np = [30, 50, 40]
L = [100., 50., 100.]
R = [3., 5., 1.5]
v = [2., 3., 4.]
@save datadir("parameters","instance_1.jld2") Np L R v

run(;param_file_name="instance_1.jld2")