using Agents
using GLMakie

#Model constants 
const repulsion_force = 1
const diameter = 1
const dt = .1
const β = 0.1
const extent = (102.,102.,102.)
const n_walkers = 1
const tMax = 5000
const dtSave = 100
const div_rate = 0.001*dt

#Create agent
mutable struct randomWalker <: AbstractAgent

    id::Int
    pos::NTuple{3,Float64}
    vel::NTuple{3,Float64}
    f::NTuple{3,Float64}

end

#Initialization
function initialize_model(n_walkers)

    space3d = ContinuousSpace(extent, extent[1])

    model = ABM(randomWalker, space3d)

    for _ in 1:n_walkers
        pos = extent./2
        vel = (0,0,0)
        f = (0,0,0)

        add_agent!(pos,model,vel,f)
    end

    return model

end

#Compute forces
function compute_forces!(walker, model)

    neighbor_ids = nearby_ids(walker,model,diameter)

    F = (0.,0.,0.)
    for id in neighbor_ids

        #Compute distance between agents
        neighbor = model[id].pos
        d = sqrt(sum((walker.pos .-neighbor).^2))

        #Add to repulsion force if closeby
        if d < diameter

            v = neighbor .- model[id].pos
            modulus = sqrt(sum(v.^2))
            if modulus > 0
                F = F .- repulsion_force .* v ./modulus
            end

        end

        walker.f = F

    end

end

#Integrator
function euler!(walker)

    # Make an euclidean step
    walker.pos = walker.pos .+ walker.vel .* dt
    walker.vel = Tuple(walker.vel .+ walker.f .*dt .- β .* walker.vel)

    #Check boundaries
    if walker.pos[1] < 0
        move_agent!(walker,(0.1,walker.pos[2],walker.pos[3]), model)
    elseif walker.pos[1] > extent[1]
        move_agent!(walker,(extent[1]-.1,walker.pos[2],walker.pos[3]), model)
    end

    if walker.pos[2] < 0
        move_agent!(walker,(walker.pos[1],0.1,walker.pos[3]), model)
    elseif walker.pos[2] > extent[2]
        move_agent!(walker,(walker.pos[1],extent[2]-.1,walker.pos[3]), model)
    end

    if walker.pos[3] < 0
        move_agent!(walker,(walker.pos[1],walker.pos[2],0.1), model)
    elseif walker.pos[3] > extent[3]
        move_agent!(walker,(walker.pos[1],walker.pos[2],extent[3]-.1), model)
    end

end

#Integrator
function division!(walker,model)

    if rand() < div_rate

        #Compute division directions
        div_direction = (randn(),randn(),randn())
        div_direction = div_direction ./ sqrt(sum(div_direction.^2))

        #Add agent
        pos_new = walker.pos .- (diameter/2)^(1/3) .* div_direction
        add_agent!(pos_new,model,walker.vel,walker.f)
        pos_new = walker.pos .+ (diameter/2)^(1/3) .* div_direction
        add_agent!(pos_new,model,walker.vel,walker.f)
        #Remove agent
        kill_agent!(walker,model)

    end

end

#Step
function agent_step!(model)

    walkers = collect(Schedulers.fastest(model))

    for id in walkers
        compute_forces!(model[id],model)
    end

    for id in walkers
        euler!(model[id])
    end
    
    for id in walkers
        division!(model[id],model)
    end

end

#Run model without saving
model = initialize_model(n_walkers)

#Plots
steps = Int64(round(tMax/dtSave))
n_steps = Int64(round(dtSave/dt))
for step in range(1,steps,step=1)

    @time step!(model,dummystep,agent_step!,n_steps)
    println("Global Time: ", step*n_steps,"/",tMax/dt)
    println("N agents: ", nagents(model),"\n")

    # x = [i[1] for i in agent[agent.step .== step,:].pos]
    # y = [i[2] for i in agent[agent.step .== step,:].pos]
    # z = [i[3] for i in agent[agent.step .== step,:].pos]
    # fig = Figure(figsize=[10,10],)
    # ax = Axis3(fig[1,1],aspect=:data)
    # meshscatter!(ax,x,y,z,markersize=diameter/2)
    # xlims!(0,extent[1])
    # ylims!(0,extent[2])
    # zlims!(0,extent[3])

    # save(string("outcome/",step,".png"),fig)

end