include("quantum_reservoir.jl")
using JSON

# Verify multithreading
println("Number of threads active is ", Threads.nthreads())

# Constants of the simulation
N_QUBITS_INPUT = 2
N_QUBITS_RESERVOIR = 5
N_QUBITS_TOTAL = N_QUBITS_INPUT + N_QUBITS_RESERVOIR
DT = 0.025
T_END = 10
N_DATAPOINTS = 10_000
OUTPUT_FOLDER = "data/experiment_run_2"


# Allowing the input states to change every time script is run
Random.seed!()
# Generating the training data
input_states = zeros(N_QUBITS_INPUT ^ 2, N_QUBITS_INPUT ^ 2, N_DATAPOINTS)im
input_state_labels = zeros(2, N_DATAPOINTS)

for ii in 1:convert(Int, N_DATAPOINTS/2)
    temp = gen_separable(N_QUBITS_INPUT ^ 2)
    input_states[:, :, ii] = temp[1]
    input_state_labels[1, ii] = temp[2] # minimum eigenvalue
    input_state_labels[2, ii] = 0 # separable = 0
end 

for ii in convert(Int, N_DATAPOINTS/2):N_DATAPOINTS
    temp = gen_entangled(N_QUBITS_INPUT ^ 2)
    input_states[:, :, ii] = temp[1]
    input_state_labels[1, ii] = temp[2] # minimum eigenvalue
    input_state_labels[2, ii] = 1 # entangled = 1
end


# Set the seed for reproducibility of reservoir configurations
Random.seed!(123)
# Reservoir configuration 
J_matrix = fully_connected_coupling_matrix(N_QUBITS_RESERVOIR)
# Random cascade terms (but fixed for each state) for symmetry breaking
gamma_C = rand(N_QUBITS_RESERVOIR)

#Initialise the reservoir
reservoir = [[0,1] for _ in 1:N_QUBITS_RESERVOIR]
interface = [1,1]
gamma_D = 0.5*sum(gamma_C)

# Timing
times1 = collect(0: DT: 2)
times2 = collect(times1[end]: DT: T_END)


# Generate data
Threads.@threads for ii in tqdm(1: N_DATAPOINTS)
    reservoir_excitations = zeros(length(times2), N_QUBITS_RESERVOIR)

    rho0 = total_rho(input_states[:, :, ii], reservoir)
    tensor_cascade = connect_interface(rho0, times1, gamma_C, gamma_D, interface)
    tensor_hamiltonian = reservoir_dynamics(tensor_cascade[:, :, end], times2, J_matrix, interface)
    
    # Measure excitations at each timestep
    for t in 1: length(times2)
        for jj in 1: N_QUBITS_RESERVOIR
            reservoir_excitations[t, jj] = abs( sig_z_exp( tensor_hamiltonian[:, :, t], (jj + N_QUBITS_INPUT) ) )
        end
    end

    # Write a dictionary to export as a JSON
    export_data = Dict(
        "reservoir_dynamics" => reservoir_excitations,
        "minimum_eigenvalue" => input_state_labels[1, ii],
        "entangled_label" => input_state_labels[2, ii]
    )

    # Write to JSON file
    open(OUTPUT_FOLDER * "/experiment_$(ii).json", "w") do f
        JSON.print(f, export_data)
    end
end
