include("quantum_reservoir.jl")

# Set the seed for reproducibility
Random.seed!(123)

# Verify multithreading
println("Number of threads active is ", Threads.nthreads())

# Defining the simulation
N_QUBITS_INPUT = 2
N_QUBITS_RESERVOIR = 5
N_QUBITS_TOTAL = N_QUBITS_INPUT + N_QUBITS_RESERVOIR
DT = 0.025
T = 10
N_DATAPOINTS = 10

# Generating the training data
input_states = zeros(N_QUBITS_INPUT ^ 2, N_QUBITS_INPUT ^ 2, N_DATAPOINTS)im
input_state_labels = zeros(2, N_DATAPOINTS)

for ii in 1:convert(Int, N_DATAPOINTS/2)
    temp = gen_separable(N_QUBITS_INPUT ^ 2)
    input_states[:, :, ii] = temp[1]
    input_state_labels[1, ii] = temp[2]
    input_state_labels[2, ii] = 0 # separable = 0
end 

for ii in convert(Int, N_DATAPOINTS/2):N_DATAPOINTS
    temp = gen_entangled(N_QUBITS_INPUT ^ 2)
    input_states[:, :, ii] = temp[1]
    input_state_labels[1, ii] = temp[2]
    input_state_labels[2, ii] = 1 # entangled = 1
end