# To call run function (see ABP output) from command line
using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")
using ArgParse

include(projectdir("src","ABP output.jl"))


#---------------------------------------------------------------------------------------------------------------------
# ARGUMENT PARSING

s = ArgParseSettings()
@add_arg_table s begin
    "--wall_condition", "-w"
        help = "wall condition for the simulation : \"periodic\", \"square\", \"ellipse\""
        arg_type = String
        default = "periodic"
    "--param_file_name", "-p" 
        help = "name of the file containing the parameters to be tested, must be stored in data/parameters !"
        arg_type = String
        required = true
    "--nb_runs", "-n"
        help = "number of runs for one parameter instance"
        arg_type = Int
        default = 1
end


#---------------------------------------------------------------------------------------------------------------------
# CALL RUNNER WITH ARGUMENTS FROM COMMAND LINE

# extract arguments
args = parse_args(s)
run(;param_file_name=args["param_file_name"], nb_runs=args["nb_runs"], wall_condition=args["wall_condition"])