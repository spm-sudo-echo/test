using DataFrames, CSV

# data from CSV
position_data = Matrix(CSV.read("/home/fantino/CELLOIDS/ABP_output/data_ellipse.csv", DataFrame))

# number of steps in the simulation
N_steps = 1000
# number of particles
N_particles = Int(size(position_data)[1] / N_steps)


# initilization for the loop
order_estimator = 0

for i_time = 0 : N_steps - 2
    # v_accumulator over particles
    v_acc_time = 0
    for i_particle = 1 : N_particles
        # the velocity is the difference between positions 
        # x (in col 3), y (in col 4) at current and next time step
        current_time = (i_time * N_particles) + i_particle
        future_time = ((i_time + 1) * N_particles) + i_particle
        println(future_time)
        vx = position_data[current_time, 3] - position_data[future_time, 3]
        vy = position_data[current_time, 4] - position_data[future_time, 4]
        # the velocity is normalized in the formula
        v_acc_time = v_acc_time .+ (vx, vy) ./ sqrt(vx^2 + vy^2)
    end 
    # accumulator over time, take the norm in the formula
    global order_estimator += sqrt(v_acc_time[1]^2 + v_acc_time[2]^2)
end

# normalize for averaging over number of particles
order_estimator /= N_particles
# normalize for averaging over number of steps
order_estimator /= N_steps

println(order_estimator)