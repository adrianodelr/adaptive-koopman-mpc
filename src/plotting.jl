# line plot colors
red_plot = RGBf(233/255, 39/255, 56/255) 
grey_plot =  RGBf(206/255, 206/255, 206/255)
yellow_plot = RGBf(235/255, 158/255, 43/255) 

function plot_tracking_results((x_akmpc,u_akmpc,t_akmpc),(x_skmpc,u_skmpc,t_skmpc),  (x_lmpc,u_lmpc,t_lmpc), (x_ref, t_ref))
    
    @assert size(x_akmpc)[1] == size(x_skmpc)[1] == size(x_lmpc)[1] == size(x_ref)[1] "Dimension of input matrices corresponding to state dimensions need to be equal."
    linew = 3

    if size(x_akmpc)[1] == 2
        fig = Figure(fontsize = 30, size=(1100, 550))
        ga = fig[1:2,1]= GridLayout() 
        
        # θ_1
        ax1b = Axis(ga[1,1], xlabel=L"t (s)",ylabel=L"$θ_1$ (rad)", width=1000,height=225,ytickformat = "{:.1f}")
        x1s = lines!(ax1b, t_skmpc[:], x_skmpc[1,:], color="gray", linewidth=linew)
        x1l = lines!(ax1b, t_lmpc[:], x_lmpc[1,:], color=yellow_plot, linewidth=linew/2)
        x1a = lines!(ax1b, t_akmpc[:], x_akmpc[1,:], color=red_plot, linewidth=linew)
        x1r = lines!(ax1b, t_ref[:], x_ref[1,:], color="black", linestyle=(:dash, :dense), linewidth=linew)
        
        # u_1 
        ax2b= Axis(ga[2,1], xlabel=L"t (s)",ylabel=L"$u_1$ (Nm)", width=1000,height=225,ytickformat = "{:.1f}")
        x1s = lines!(ax2b, t_skmpc[1:end-1], u_skmpc[1,:], color="gray", linewidth=linew)
        x1l = lines!(ax2b, t_lmpc[1:end-1], u_lmpc[1,:], color=yellow_plot, label=L"$$linearization",linewidth=linew, )
        x1a = lines!(ax2b, t_akmpc[1:end-1], u_akmpc[1,:], color=red_plot, linewidth=linew)

        hidexdecorations!(ax1b, ticks = true, grid=false)
        rowgap!(ga, 0)
        colgap!(ga, 0)

        axislegend(ax1b, [x1r, x1l,x1a, x1s], [L"$$reference", L"$$linearization MPC",L"$$adaptive KMPC", L"$$static KMPC"], position = :rb)

    elseif size(x_akmpc)[1] == 4
        fig = Figure(fontsize = 30, size=(1100, 1000))
        ga = fig[1:4,1]= GridLayout()

        # θ_1
        ax1b = Axis(ga[1,1], xlabel=L"t (s)",ylabel=L"$θ_1$ (rad)", width=1000,height=225,ytickformat = "{:.1f}")
        x1l = lines!(ax1b, t_lmpc[:], x_lmpc[1,:], color=yellow_plot, linewidth=linew)
        x1s = lines!(ax1b, t_skmpc[:], x_skmpc[1,:], color="gray", linewidth=linew)
        x1a = lines!(ax1b, t_akmpc[:], x_akmpc[1,:], color=red_plot, linewidth=linew)
        x1r = lines!(ax1b, t_ref[:], x_ref[1,:], color="black", linestyle=(:dash, :dense), linewidth=linew)
        
        # θ_2
        ax1a = Axis(ga[2,1], xlabel=L"$t$ (s)",ylabel=L"$θ_2$ (rad)", width=1000,height=225,ytickformat = "{:.1f}")
        x1l = lines!(ax1a, t_lmpc[:], x_lmpc[2,:], color=yellow_plot,linewidth=linew)
        x1s = lines!(ax1a, t_skmpc[:], x_skmpc[2,:], color="gray", linewidth=linew)
        x1a = lines!(ax1a, t_akmpc[:], x_akmpc[2,:], color=red_plot, linewidth=linew)
        x1r = lines!(ax1a, t_ref[:], x_ref[2,:], color="black", linestyle=(:dash, :dense), linewidth=linew)
        
        # u_1 
        ax2b= Axis(ga[3,1], xlabel=L"t (s)",ylabel=L"$u_1$ (Nm)", width=1000,height=225,ytickformat = "{:.1f}")
        x1a = lines!(ax2b, t_akmpc[1:end-1], u_akmpc[1,:], color=red_plot,linewidth=linew)
        x1s = lines!(ax2b, t_skmpc[1:end-1], u_skmpc[1,:], color="gray",linewidth=linew)
        x1l = lines!(ax2b, t_lmpc[1:end-1], u_lmpc[1,:], color=yellow_plot,linewidth=linew)
        
        # u_2 
        ax2a= Axis(ga[4,1], xlabel=L"$t$ (s)",ylabel=L"$u_2$ (Nm)", width=1000,height=225,ytickformat = "{:.1f}")
        x1a = lines!(ax2a, t_akmpc[1:end-1], u_akmpc[2,:], color=red_plot, linewidth=linew)
        x1s = lines!(ax2a, t_skmpc[1:end-1], u_skmpc[2,:], color="gray",linewidth=linew)
        x1l = lines!(ax2a, t_lmpc[1:end-1], u_lmpc[2,:], color=yellow_plot,linewidth=linew)
        
        axislegend(ax1b, [x1r, x1l,x1a, x1s], [L"$$reference", L"$$linearization MPC",L"$$adaptive KMPC", L"$$static KMPC"], position = :rb)
        
        hidexdecorations!(ax1a, ticks = true, grid=false)
        hidexdecorations!(ax1b, ticks = true, grid=false)
        hidexdecorations!(ax2b, ticks = true, grid=false)


        rowgap!(ga, 0)
    else    
        println("Dimension of input matrices corresponding to state dimensions can either be two or four.")
    end 

    resize_to_layout!(fig)
    return fig
end  



function plot_performance_metrics((x_akmpc,u_akmpc),(x_skmpc,u_skmpc), (x_lmpc,u_lmpc), x_ref, Δt )

    @assert size(x_akmpc)[1] == size(x_skmpc)[1] == size(x_lmpc)[1] == size(x_ref)[1] "Dimension of input matrices corresponding to state dimensions need to be equal."
    
    m_dof = 0.5size(x_akmpc,1) |> Int 
    N = size(x_ref,2)

    lmpc_power_pos_total,  lmpc_power_neg_total  = 0, 0
    akmpc_power_pos_total, akmpc_power_neg_total = 0, 0
    skmpc_power_pos_total, skmpc_power_neg_total = 0, 0

    lmpc_MSE  = zeros(m_dof)
    akmpc_MSE = zeros(m_dof)
    skmpc_MSE = zeros(m_dof)
    
    for i in 1:m_dof 
        # power computation 
        # lmpc 
        lmpc_power_pos =  max.(0,u_lmpc[i,:].*x_lmpc[i+m_dof,1:end-1])                                                       
        lmpc_power_neg =  min.(0,u_lmpc[i,:].*x_lmpc[i+m_dof,1:end-1])                                                      

        lmpc_power_pos_total += 0.5*(lmpc_power_pos[2:end]+lmpc_power_pos[1:end-1]).*Δt |> sum |> abs             
        lmpc_power_neg_total += 0.5*(lmpc_power_neg[2:end]+lmpc_power_neg[1:end-1]).*Δt |> sum |> abs             
        # akmpc 
        akmpc_power_pos =  max.(0,u_akmpc[i,:].*x_akmpc[i+m_dof,1:end-1])                                                       
        akmpc_power_neg =  min.(0,u_akmpc[i,:].*x_akmpc[i+m_dof,1:end-1])                                                        
        akmpc_power_pos_total += 0.5*(akmpc_power_pos[2:end]+akmpc_power_pos[1:end-1])*Δt |> sum |> abs             
        akmpc_power_neg_total += 0.5*(akmpc_power_neg[2:end]+akmpc_power_neg[1:end-1])*Δt |> sum |> abs              
        # skmpc 
        skmpc_power_pos =  max.(0,u_skmpc[i,:].*x_skmpc[i+m_dof,1:end-1])                                                       
        skmpc_power_neg =  min.(0,u_skmpc[i,:].*x_skmpc[i+m_dof,1:end-1])                                                      
        skmpc_power_pos_total += 0.5*(skmpc_power_pos[2:end]+skmpc_power_pos[1:end-1])*Δt |> sum |> abs             
        skmpc_power_neg_total += 0.5*(skmpc_power_neg[2:end]+skmpc_power_neg[1:end-1])*Δt |> sum |> abs            
        
        # mean squared error
        lmpc_MSE[i]  = 1/N*sum((x_lmpc[i,:]-x_ref[i,:]).^2)
        akmpc_MSE[i] = 1/N*sum((x_akmpc[i,:]-x_ref[i,:]).^2)         
        skmpc_MSE[i] = 1/N*sum((x_skmpc[i,:]-x_ref[i,:]).^2)          
    end 

    colors = [yellow_plot, red_plot, "gray"]
    positions = [1,1,1, 2,2,2]
    dodge = [1,2,3, 1,2,3]

    if m_dof == 1
        fig = Figure(fontsize = 30, size=(800, 650))
        ga = fig[1,1:2] = GridLayout()
    
        ax1a = Axis(ga[1,1], width = 400, height = 500, xticks = (1:2, [L"$$Pos.", L"$$Neg."]), title = L"$$Power (J)")
        height = [lmpc_power_pos_total,
                akmpc_power_pos_total, 
                skmpc_power_pos_total, 
                lmpc_power_neg_total,
                akmpc_power_neg_total, 
                skmpc_power_neg_total]

        barplot!(ax1a, positions, height,
                dodge = dodge,
                color = colors[dodge])


        ax1b  = Axis(ga[1,2], width = 175, height = 500, xticks = (1:1, [L"$θ_1$"]), yaxisposition = :right, title = L"$$MSE")
        height = [lmpc_MSE[1], akmpc_MSE[1], skmpc_MSE[1]]

        barplot!(ax1b, positions[1:3], height,
                dodge = dodge[1:3],
                color = colors[dodge[1:3]])

        colgap!(ga, 0)

    elseif m_dof == 2
        fig = Figure(fontsize = 30, size=(1000, 650))
        ga = fig[1,1:2] = GridLayout()

        ax1a = Axis(ga[1,1], width = 400, height = 500, xticks = (1:2, [L"$$Pos.", L"$$Neg."]), title = L"$$Power (J)")
        height = [lmpc_power_pos_total,
                akmpc_power_pos_total, 
                skmpc_power_pos_total, 
                lmpc_power_neg_total,
                akmpc_power_neg_total, 
                skmpc_power_neg_total]

        barplot!(ax1a, positions, height,
                dodge = dodge,
                color = colors[dodge])

        ax1b  = Axis(ga[1,2], width = 400, height = 500, xticks = (1:2, [L"$θ_1$", L"$θ_2$"]), yaxisposition = :right, title = L"$$MSE")

        height = [lmpc_MSE[1], akmpc_MSE[1], skmpc_MSE[1], lmpc_MSE[2], akmpc_MSE[2], skmpc_MSE[2]]

        barplot!(ax1b, positions, height,
                dodge = dodge,
                color = colors[dodge])
        
        colgap!(ga, 0)
    end 
    # Legend
    labels = [L"$$lin. MPC", L"$$adapt. KMPC", L"$$stat. KMPC"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    title = L"$$Controller"    
    Legend(fig[2,1:2], elements, labels, title,  orientation = :horizontal, titleposition = :left)


    fig 

    return fig
end 