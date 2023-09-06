# calculate the velocity and velocity polarization of particles , i.e., velocity order parameter
# Date: 29-05-2023
# Method: Absoulte velocity calculation at each instant
# Input: File from the output.jl, which has positions and orientation at each instant
# Output: Velocity of each particle at every instant and average velocity of the ensemble

using  Plots, LaTeXStrings, Statistics, CSV, DataFrames,CategoricalArrays, FilePathsBase
import PackageName


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Estimation of the velocity based on the simulation output file (must be in data/sims !)
# 
function estimate_velocity(source_file::String, averaging::Bool=false, save::Bool=false)
    df= CSV.read(joinpath(source_file), DataFrame)
    # steps= df[!, :N]
    # time= df[!,:Time]
    # x= df[!,:xpos]
    # y= df[!,:ypos]
    df[!,:N] = categorical(df[!,:N],compress=true)  # it sorts out time step data 
    ## Group dataframe by values in categorical column
    gdf = groupby(df,:N,sort=true)                  # only 1000 data groups because I have omitted 100 time steps means 1 s
    # averaging or not
    if(!averaging)
        vel_x = [diff(g[!, :xpos]) for g in gdf]    # x velocity [particle i] [time step j]
        vel_y = [diff(g[!,:ypos]) for g in gdf]     # same for y
    else
        vel_x = [mean(diff(g[!, :xpos])) for g in gdf]   
        vel_y = [mean(diff(g[!,:ypos])) for g in gdf]     
    end
    # exportation or not
    if(!save)
        return vel_x, vel_y
    else
        velocity_file_name = filename(AbstractPath(source_file)) 
        # the name will be different if we average or not
        averaging ? velocity_file_name *= "-mean_velocity.csv" : velocity_file_name *= "-velocity.csv" 
        # velocity_path
        # open(data_vel_path, "w")
        # data = DataFrame(
        #     N= pnumber,
        #     Time= time,
        #     xpos= x,
        #     ypos= y,
        #     orientation=θ) 
        # CSV.write(data_vel_path, data)
    end
end


function polarization_factor(source_file_path::String)
    # vel_magnitude = [sqrt.(vel_x[i].^2 + vel_y[i].^2) for i in eachindex(vel_x)]  # velocity magnitude at each time second
    # vp = (vel_x./vel).+ (vel_y./vel)            # polarization factor 

end

#plot(vel[100],vel[1000])
