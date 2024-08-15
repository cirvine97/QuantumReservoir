using LinearAlgebra
using QuantumInformation
using Plots
using ProgressBars
using DataFrames
using CSV
using Random

# Set the seed for reproducibility
Random.seed!(123)


function fully_connected_coupling_matrix(N)
    """Return the coupling matrix for a fully connected reservoir where all N qubits interact.
    
        Args:
            N (int): Number of qubits in reservoir
        Returns:
            J (Matrix): Coupling matrix
    """
    J = rand(N, N)

    # Set diagonal elements for zero (no self-interaction)
    for i in 1:N
        J[i, i] = 0
    end

    return J
end


function commutator(A,B,kind)
    comm=0
    if kind=="anti"
        comm=A*B+B*A
    elseif kind=="normal"
        comm=A*B-B*A
    end
end


function total_rho(interface_rho, reservoir)
    """Return the initial density matrix of the system given the interface as a dnsity matrix"""
    placeholder=reservoir[1]
    for ii in 2:length(reservoir)
        placeholder=kron(placeholder, reservoir[ii])
    end
    rhoR=placeholder*placeholder' #' is the conjugate transpose in Julia
    #Final tensor 
    rho=kron(interface_rho, rhoR)
end 


function ptb(rho)
    """Return the partial transpose of rho about subsystem B for a 4x4 density matrix"""
    a=zeros(Complex{Float64},(length(rho[:,1]), length(rho[:,1])))
    a[1:2,1:2]=transpose(rho[1:2,1:2])
    a[1:2,3:4]=transpose(rho[1:2,3:4])
    a[3:4,1:2]=transpose(rho[3:4,1:2])
    a[3:4,3:4]=transpose(rho[3:4,3:4])
    return a
end


function tr2(rho)
    """Return the purity of rho"""
    pur=real(tr(rho^2))
end


function linear_entropy(rho)
    """Calculate the linear entropy of the given density matrix"""
    N=length(rho[:,1])
    le=(N/(N-1)) * (1-tr2(rho))
end


function gen_purestate(N)
    """Return an NxN density matrix for a random pure state"""
    #Generatte a random complex Nx1 vector 
    p=rand(ComplexF64, N) - rand(ComplexF64, N)
    #Normalise the vector 
    normalize!(p)
    #Generate the density matrix
    dens=p*p'
end


function gen_mixedstate(N)
    """Return an NxN density matrix for a random mixed state"""
    m=rand(Complex{Float64}, (N,N))-rand(Complex{Float64}, (N,N))
    #Normalise 
    M_dens=(m*m')/tr(m*m')
end 


function gen_entangled(N)
    """Return an NxN density matrix for a random entangled system"""
    #Define local variables to function
    temp=Array{ComplexF64,(N,N)}[]
    eigenvalues=zeros(N)
    
    while minimum(eigenvalues)>=0
        #Generate a random mixedstate
        temp=gen_mixedstate(N)
        #Perform the partial transpose 
        ptm=ptb(temp)
        #Find the eigenvalues 
        eigenvalues=real(eigvals(ptm))
    end 
    return temp, minimum(eigenvalues), linear_entropy(temp)
end


function gen_pure_entangled(N)
    """Generate a pure NxN entangled state"""
    temp=zeros(N,N)im
    eigenvalues=zeros(N)
    #Enforce the eigenvalues of the PT must be negative so system is entangled 
    while minimum(eigenvalues)>=0
        #Generate the random purestate
        temp=gen_purestate(N)
        PT=ptb(temp)
        eigenvalues=real(eigvals(PT))
    end
    return temp, minimum(eigenvalues), linear_entropy(temp)
end


function gen_separable(N)
    """Return a separable density matrix"""
    #Define the local variables 
    temp=Array{ComplexF64, (N,N)}[]
    eigenvalues=zeros(N)
    
    #Generate the random mixedstate
    temp=gen_mixedstate(N)
    #Find the partial transpose 
    ptm=ptb(temp)
    #Find the eigenvalues of the PT
    eigenvalues=real(eigvals(ptm))
    #Enforce that the eigenvalues of the PT must be positive
    while minimum(eigenvalues)<=0
        p=rand(Float64) #Random number (0,1]
        temp=p*temp + 0.25*(1-p)*I(4)
        ptm=ptb(temp)
        eigenvalues=real(eigvals(ptm))
    end
    return temp, minimum(eigenvalues), linear_entropy(temp)
end


function gen_pure_separable(N)
    """Return a separable density matrix beginning from a pure state"""
    #Define the local variables 
    temp=Array{ComplexF64, (N,N)}[]
    eigenvalues=zeros(N)
    
    #Generate the random mixedstate
    temp=gen_purestate(N)
    #Find the partial transpose 
    ptm=ptb(temp)
    #Find the eigenvalues of the PT
    eigenvalues=real(eigvals(ptm))
    #Enforce that the eigenvalues of the PT must be positive
    while minimum(eigenvalues)<=0
        p=rand(Float64) #Random number (0,1]
        temp=p*temp + 0.25*(1-p)*I(4)
        ptm=ptb(temp)
        eigenvalues=real(eigvals(ptm))
    end
    return temp, minimum(eigenvalues), linear_entropy(temp)
end 
    
    
function cascade(rho, interface_vector, gamma_vector)
    """Generate the cascade term for the dynamics from interface to specified reservoir qubits"""
    #Define the needed Pauli Matrices first 
    global sigmam=[0 0 ; 1 0]
    global sigmap=[0 1 ; 0 0]
    #Define the local variables
    interface_m=zeros(2,2)
    interface_p=zeros(2,2)
    N_int=length(interface_vector)
    sig_I_minus=Matrix{ComplexF64}
    sig_I_plus=Matrix{ComplexF64}
    sig_pos_list=Matrix{ComplexF64}
    sig_minus_list=Matrix{ComplexF64}
    cas=zeros(length(rho[:,1]),length(rho[:,1]))
    
    #Interface vector is boolean indicating if interface qubit should connect to reservoir
    if interface_vector[1]==0
        interface_m=I(2) #Identity matrix
        interface_p=I(2)
    elseif interface_vector[1]==1
        interface_m=sigmam
        interface_p=sigmap
    end 
    #Populate the rest of the matrices for the operators 
    for ii in 2:length(interface_vector)
        if interface_vector[ii]==0
            interface_m=interface_m⊗I(2)
            interface_p=interface_p⊗I(2)
        elseif interface_vector[ii]==1
            interface_m=interface_m⊗sigmam
            interface_p=interface_p⊗sigmap
        end
    end
    
    #Begin on the reservoir terms 
    #Include a clause for the 1 qubit case 
    if length(gamma_vector)==1
        sig_I_minus=interface_m⊗I(2)
        sig_I_plus=interface_p⊗I(2)
        sig_1_minus=I(2^N_int)⊗sigmam
        sig_1_plus=I(2^N_int)⊗sigmap
        cas=gamma_vector[1]*(commutator(sig_I_minus*rho,sig_1_plus,"normal")+
        commutator(sig_1_minus,rho*sig_1_plus,"normal"))
    
    else
        #Generate the needed pauli operators in system size 
        sig_I_minus=interface_m⊗I(convert(Int,(length(rho[:,1])/(2^N_int))))
        sig_I_plus=interface_p⊗I(convert(Int,(length(rho[:,1])/(2^N_int))))
        sig_pos_list=zeros(2^(length(gamma_vector)+N_int),
        2^(length(gamma_vector)+N_int),
        length(gamma_vector))
        sig_minus_list=zeros(2^(length(gamma_vector)+N_int),
        2^(length(gamma_vector)+N_int),
        length(gamma_vector))
        
        #Populate the operator list for the reservoir with the first and last elements first 
        sig_pos_list[:,:,1]=I(2^N_int)⊗(sigmap⊗I(2^(length(gamma_vector)-1)))
        sig_minus_list[:,:,1]=I(2^N_int)⊗(sigmam⊗I(2^(length(gamma_vector)-1)))
        sig_pos_list[:,:,end]=I(2^N_int)⊗(I(2^(length(gamma_vector)-1))⊗sigmap)
        sig_minus_list[:,:,end]=I(2^N_int)⊗(I(2^(length(gamma_vector)-1))⊗sigmam)
        #Populate intermediary terms if needed 
        if length(gamma_vector)>2
            for ii in 2:(length(gamma_vector)-1)
                sig_pos_list[:,:,ii]=I(2^N_int)⊗I(2^(ii-1))⊗sigmap⊗I(2^(length(gamma_vector)-ii))
                sig_minus_list[:,:,ii]=I(2^N_int)⊗I(2^(ii-1))⊗sigmam⊗I(2^(length(gamma_vector)-ii))
            end
        end

        #Populate the cascade terms for this situation
        for ii in 1:length(gamma_vector)            
            cas+=gamma_vector[ii]*(commutator(sig_I_minus*rho, sig_pos_list[:,:,ii], "normal")+
            commutator(sig_minus_list[:,:,ii], rho*sig_I_plus, "normal"))
        end
    end
    return cas
end
        

function dissipate(rho, interface_vector, gamma_D)
    """Dissipation for each interface qubit to ensure physicality"""
    #Define the local variables 
    interface_m=Matrix{ComplexF64}
    interface_p=Matrix{ComplexF64} 
    sig_I_minus=Matrix{ComplexF64} 
    sig_I_plus=Matrix{ComplexF64} 
    diss=zeros(length(rho[:,1]), length(rho[:,1]))
    N_int=length(interface_vector)
    
    #Generate the interface operators 
    if interface_vector[1]==0
        interface_m=I(2) #Identity matrix
        interface_p=I(2)
    elseif interface_vector[1]==1
        interface_m=sigmam
        interface_p=sigmap
    end 
    #Populate the rest of the matrices for the operators 
    for ii in 2:length(interface_vector)
        if interface_vector[ii]==0
            interface_m=interface_m⊗I(2)
            interface_p=interface_p⊗I(2)
        elseif interface_vector[ii]==1
            interface_m=interface_m⊗sigmam
            interface_p=interface_p⊗sigmap
        end
    end
    
    #Define the needed Pauli Operators
    sig_I_minus=interface_m⊗I(convert(Int, (length(rho[:,1])/(2^N_int))))
    sig_I_plus=interface_p⊗I(convert(Int, (length(rho[:,1])/(2^N_int))))
    
    #Generate the dissipation term using anticommutators for the Lindblad equation
    diss=gamma_D*(2*sig_I_minus*rho*sig_I_plus - commutator(sig_I_plus*sig_I_minus, rho, "anti"))
end


function hamiltonian_R(rho, J_matrix)
    """Return the Hamiltonian of the reservoir given the coupling between qubits"""
    #Initialise the local variables
    N=convert(Int,log2(length(rho[:,1]))) #Dimension of the whole system
    N_res=convert(Int,length(J_matrix[:,1])) #Dimension of the reservoir 
    N_int=convert(Int,log2(length(rho[:,1]))-length(J_matrix[:,1])) #Dimension of the interface
    H_R=zeros(2^N_res,2^N_res) 
    sig_pos_list=zeros(2^N_res, 2^N_res, N_res) #Operators for the Hamiltonian
    sig_minus_list=zeros(2^N_res, 2^N_res, N_res)
    
    #Populate the operator list 
    sig_pos_list[:,:,1]=sigmap⊗I(2^(N_res-1))
    sig_minus_list[:,:,1]=sigmam⊗I(2^(N_res-1))
    sig_pos_list[:,:,end]=I(2^(N_res-1))⊗sigmap
    sig_minus_list[:,:,end]=I(2^(N_res-1))⊗sigmam
    #Any interior terms
    if N_res>2
        for ii in 2:(N_res-1)
            sig_pos_list[:,:,ii]=I(2^(ii-1))⊗sigmap⊗I(2^(N_res-ii))
            sig_minus_list[:,:,ii]=I(2^(ii-1))⊗sigmam⊗I(2^(N_res-ii))
        end
    end
    
    #Populate the hamiltonian matrix
    for ii in 1:N_res
        for jj in 1:N_res
            H_R+=J_matrix[ii,jj]*(sig_pos_list[:,:,ii]*sig_minus_list[:,:,jj]+sig_pos_list[:,:,jj]*sig_minus_list[:,:,ii])
        end
    end
    
    #Now tensor with identity of interface to match system size
    H_R=I(2^N_int)⊗H_R
end


function connect_interface(rho0, times, gamma_C, gamma_D, interface_vector)
    """Connect the interface and reservoir and return density tensor in time"""
    #Initialise the local variables 
    dens_tensor=zeros(length(rho0[:,1]), length(rho0[:,1]), length(times))im
    #Set the first term 
    dens_tensor[:,:,1]=rho0
    #Find the time separation 
    dt=times[2]-times[1]

    #Evolve using Runge-Kutta O(dt^4) method
    for ii in 2:length(times)
        #Determine the Runge-Kutta terms 
        k1=(cascade(dens_tensor[:,:,ii-1], interface_vector, gamma_C)+
            dissipate(dens_tensor[:,:,ii-1], interface_vector, gamma_D))*dt
        k2=(cascade((dens_tensor[:,:,ii-1]+0.5*k1), interface_vector, gamma_C)+
            dissipate((dens_tensor[:,:,ii-1]+0.5*k1), interface_vector, gamma_D))*dt
        k3=(cascade((dens_tensor[:,:,ii-1]+0.5*k2), interface_vector, gamma_C)+
            dissipate((dens_tensor[:,:,ii-1]+0.5*k2), interface_vector, gamma_D))*dt
        k4=(cascade((dens_tensor[:,:,ii-1]+k3), interface_vector, gamma_C)+
            dissipate((dens_tensor[:,:,ii-1]+k3), interface_vector, gamma_D))*dt
        #Solve the Taylor expansion
        dens_tensor[:,:,ii]=dens_tensor[:,:,ii-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    end 
    return dens_tensor
end


function reservoir_dynamics(rho0, times, J_matrix, interface_vector)
    """Evolve the dynamics inside the reservoir after transfer of excitation"""
    #Initialise the local variables 
    dens_tensor=zeros(length(rho0[:,1]), length(rho0[:,1]), length(times))im
    #Set the first term 
    dens_tensor[:,:,1]=rho0
    #Find the time separation 
    dt=times[2]-times[1]
    
    #Solve the SE with Runge-Kutta O(dt^4) methods
    for ii in 2:length(times)
        #Determine the Runge-Kutta terms 
        k1=(-im * commutator(hamiltonian_R(dens_tensor[:,:,ii-1], J_matrix), dens_tensor[:,:,ii-1], "normal"))*dt
        k2=(-im * commutator(hamiltonian_R(dens_tensor[:,:,ii-1], J_matrix), (dens_tensor[:,:,ii-1]+0.5*k1), "normal"))*dt
        k3=(-im * commutator(hamiltonian_R(dens_tensor[:,:,ii-1], J_matrix), (dens_tensor[:,:,ii-1]+0.5*k2), "normal"))*dt
        k4=(-im * commutator(hamiltonian_R(dens_tensor[:,:,ii-1], J_matrix), (dens_tensor[:,:,ii-1]+k3), "normal"))*dt
        #Solve the Taylor expansion 
        dens_tensor[:,:,ii]=dens_tensor[:,:,ii-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    end
    return dens_tensor
end


function reduced_density_matrix(rho, qubit)
    """Return the reduced density matrix for the wanted qubit out of rho"""
    #Total number of qubits of the system 
    N=length(rho[:,1])
    #Find the dimension of the system leading up to the wanted qubit
    N_before=2^(qubit-1)
    #The dimension of the system after the wanted qubit 
    N_after=convert(Int, N/2^(qubit))
    #Find the reduced density matrix
    rd=ptrace(rho, [N_before,2,N_after],[1,3])
end


function sig_z_exp(rho, qubit)
    """Give the expectation value of σz for a given qubit"""
    #Define the needed Pauli Operator 
    sigmaz=[1 0 ;
            0 -1]
    #Find the reduced density matrix of the wanted qubit 
    rho_q=reduced_density_matrix(rho, qubit)
    #Find the expectation value 
    exp=real(tr(sigmaz*rho_q))
end


function reservoir_excitations(N, J_matrix, gamma_C, gamma_D, cascade_time, reservoir_time, dt, reservoir_state, interface)
    global data=zeros(4,4,N)im
    global labels=zeros(2,N)

    #Eigenvalue and linear entropy will both be computed as labels 
    for ii in 1:convert(Int,N/4)
        temp=gen_separable(4)
        data[:,:,ii]=temp[1]
        labels[1,ii]=temp[2]
        labels[2,ii]=temp[3]
    end
    for ii in convert(Int, N/4):convert(Int,N/2)
        temp=gen_pure_separable(4)
        data[:,:,ii]=temp[1]
        labels[1,ii]=temp[2]
        labels[2,ii]=temp[3]
    end
    for ii in convert(Int,N/2):convert(Int,3N/4)
        temp=gen_entangled(4)
        data[:,:,ii]=temp[1]
        labels[1,ii]=temp[2]
        labels[2,ii]=temp[3]
    end
    for ii in convert(Int,3N/4):N
        temp=gen_pure_entangled(4)
        data[:,:,ii]=temp[1]
        labels[1,ii]=temp[2]
        labels[2,ii]=temp[3]
    end
    
    #Timing
    times1=collect(0:dt:cascade_time)
    times2=collect(times1[end]:dt:reservoir_time)
    times=cat(times1,times2,dims=1)
    
    #Tensor to hold the excitations of each qubit for each interation of the reservoir
    excitations=zeros(N, length(reservoir_state))
    
    #Generate
    Threads.@threads for ii in tqdm(1:N)
        #Find the initial density matrix 
        rho0=total_rho(data[:,:,ii], reservoir_state)
        #find the cascade term 
        tensor_cascade=connect_interface(rho0, times1, gamma_C, gamma_D, interface) 
        #Allow reservoir dynamics after 
        tensor_hamiltonian=reservoir_dynamics(tensor_cascade[:,:,end], times2, J_matrix, interface) 
        #Measure the excitations each time 
        for jj in 1:length(reservoir_state)
            excitations[ii,jj]=sig_z_exp(tensor_hamiltonian[:,:,end], (jj+2))
        end
    end

    return excitations
end