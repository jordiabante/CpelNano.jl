using CpelNano
using Test

@testset "Computations check" begin

    ## Marginal expectations 
    
    # Transfer matrix objects
    @test CpelNano.get_u(0.0)==[1.0,1.0]
    @test CpelNano.get_W(0.0,0.0,0.0)==[1.0 1.0;1.0 1.0]

    # Partition function
    L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    u1 = CpelNano.get_u(αs[1]); uN = CpelNano.get_u(αs[end]);
    Ws = [CpelNano.get_W(αs[n],αs[n+1],βs[n]) for n=1:length(βs)];
    Z = CpelNano.get_Z(u1,uN,Ws)

    # Partition function
    @test Z==1024.0
    
    # Expected value of X
    @test CpelNano.get_E_X(u1,uN,Ws,Z)==zeros(10)
    
    # Expected value of XX
    @test CpelNano.get_E_XX(u1,uN,Ws,Z)==zeros(9)
    
    # Log-sum trick
    L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    logu1 = CpelNano.get_log_u(αs[1]); 
    loguN = CpelNano.get_log_u(αs[end]);
    logWs = [CpelNano.get_log_W(αs[n],αs[n+1],βs[n]) for n=1:length(βs)];
    
    # Partition function
    logZ = CpelNano.get_log_Z(logu1,loguN,logWs)
    @test logZ==log(1024)
    
    # Expected value of X
    @test isapprox(CpelNano.get_E_X_log(logu1,loguN,logWs,logZ),zeros(10);atol=1e-5)
    
    # Expected value of XX
    @test isapprox(CpelNano.get_E_XX_log(logu1,loguN,logWs,logZ),zeros(9);atol=1e-5)
    
    ## Conditional expectations 
    
    # Regular transfer matrix computation
    L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    obs = fill(CpelNano.MethCallCpgGrp(-50.0,-50.0),L);
    u1 = CpelNano.get_uc(αs[1],obs[1]); uN = CpelNano.get_uc(αs[end],obs[end]);
    Ws = [CpelNano.get_Wc(αs[n],αs[n+1],βs[n],obs[n],obs[n+1]) for n=1:length(βs)];
    
    # Partition function
    Zc = CpelNano.get_Zc(u1,uN,Ws);
    @test Zc==7.295566240503076e-215
    
    # Expected value of X
    @test CpelNano.get_Ec_X(u1,uN,Ws,Zc)==zeros(10)
    
    # Expected value of XX
    @test CpelNano.get_Ec_XX(u1,uN,Ws,Zc)==zeros(9)
    
    # Log-sum trick
    L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    obs = fill(CpelNano.MethCallCpgGrp(-100.0,-100.0),L);
    logu1 = CpelNano.get_log_uc(αs[1],obs[1]); 
    loguN = CpelNano.get_log_uc(αs[end],obs[end]);
    logWs = [CpelNano.get_log_Wc(αs[n],αs[n+1],βs[n],obs[n],obs[n+1]) for n=1:length(βs)];
    
    # Partition function
    logZc = CpelNano.get_log_Zc(logu1,loguN,logWs)
    @test logZc==-993.0685281944009
    
    # Expected value of X
    @test isapprox(CpelNano.get_Ec_X_log(logu1,loguN,logWs,logZc),zeros(10);atol=1e-5)
    
    # Expected value of XX
    @test isapprox(CpelNano.get_Ec_XX_log(logu1,loguN,logWs,logZc),zeros(9);atol=1e-5)

end
