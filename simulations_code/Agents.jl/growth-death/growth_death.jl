using Agents
using GLMakie

#Model constants
const repulsion_force = 1
const diameter = 1
const dt = .1
const β = 0.1
const extent = (72.,72.,72.)
const extent_out = (102.,102.,102.)
const n_walkers = 15
const n_steps = 100
const save_step = 1
const div_rate = 0.1*dt
const death_rate = 0.08*dt
const plot = false

#Create agent
mutable struct randomWalker <: AbstractAgent
    id::Int
    pos::NTuple{3,Float64}
    vel::NTuple{3,Float64}
    f::NTuple{3,Float64}
end

#Initialization
function initialize_model(n_walkers)
    space2d = ContinuousSpace(extent_out, 3*diameter)

    model = ABM(randomWalker, space2d)

    for _ in 1:n_walkers
        pos = (extent[1]*rand()+extent_out[1]/4., 
                extent[2]*rand()+extent_out[2]/4., 
                extent[3]*rand()+extent_out[3]/4.)
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
        d = sqrt(sum((walker.pos .- neighbor).^2))

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

            pos_new = walker.pos .- diameter/sqrt(3) .* div_direction
            move_agent!(walker, walker.pos .+ diameter/sqrt(3) .* div_direction, model)

            #Add agent
            add_agent!(pos_new,model,walker.vel,walker.f)
            
        elseif rand() < death_rate

            #Add agent
            kill_agent!(walker,model)

        end

end

#step
function agent_step!(model)
    for id in Schedulers.fastest(model)
        compute_forces!(model[id],model)
    end
    for id in Schedulers.fastest(model)
        euler!(model[id])
    end
    for id in Schedulers.fastest(model)
        division!(model[id],model)
    end
end

#Run model without saving
model = initialize_model(n_walkers)
#run!(model,dummystep,step!,n_steps,adata=[:pos],when=range(1,n_steps,step=save_step))

steps = Int64(round(n_steps/save_step))
for step in range(1,steps,step=1)
    @time step!(model,dummystep,agent_step!,save_step)

    #Plots
    if plot
        x = [i[1] for i in [agent.pos for agent in allagents(model)]]
        y = [i[2] for i in [agent.pos for agent in allagents(model)]]
        z = [i[3] for i in [agent.pos for agent in allagents(model)]]
        fig = Figure(figsize=[10,10],)
        ax = Axis3(fig[1,1],aspect=:data)
        meshscatter!(ax,x,y,z,markersize=diameter/2)
        xlims!(0,extent[1])
        ylims!(0,extent[2])
        zlims!(0,extent[3])

        save(string("outcome/",step,".png"),fig)
    end
end