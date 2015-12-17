# 15-462-Final-Project
The final project for 15-462

Zach Shearer (zshearer) and Chris Kaffine (ckaffine) (partner project)

##Proposal

Our goal is to be able to make high quality, efficient renderings of real-world physical simulations  with PDEs on meshes, by combining raytracing with mesh dynamics. Specifically, we'd like to be able to make realistic surface tensioning renderings.

1. We will speed up bvh construction with parallelism, as specified in option F.
2. We'll implement bidirectional pathtracing.
3. We'll implement surface tensioning based on the methods described in the CGL paper linked in the writeup.
4. We'll implement a few extra liquid materials, such as water.
5. Add some functionality to the GUI to allow the user to enter disturbances to the system

##Implementation

1. We were able to successfully implement improved BVH construction. Below are some statistics showing the performance improvement with the new algorithm
2. We attempted to implement bidirectional pathtracing, but were unable to finish in time. Below are some images created with the most recent version
3. The surface tension algorithm we attempted to implement consists of two steps: a) Perform a volume preserving mean curvature flow on the surface, then b) use the distance between this smoother surface and the original surface as the initial displacement in a wave equation simulation. 
	a. This step turned out to be significantly more challenging than anticipated, so we did not end up using it for our main results. On any mesh more intricate than the quad ball, the curvature flow ends up severely degrading the quality of mesh after a reasonable number of iterations. We attempted to alleviate this by implementing a remeshing scheme (based on [this paper] (https://www.cs.sfu.ca/~haoz/pubs/dyer_et_al_sgp07.pdf)) to ensure that the mesh remains Delaunay throughout the simulation. This, combined with a mesh reduction scheme that collapses sufficiently small edges, was able to improve the results, but not enough to be usable in our final results. However, there are some example videos of our implementation in the handin folder.
	b. We successfully implemented the wave equation dynamics, and were able to use this to generate our final results.
4. We successfully implemented a BSDF for glossy materials, which we used to create a BSDF for water.
5. We ended up adding some rudimentary functionality to the GUI to enable the set up and rendering of liquid animations. The extensions are admittedly not integrated perfectly, but the functionality is summarized below. Also, while in animation mode, you can click and drag vertices along their normals to generate displacements to be used in the wave equation simulation.

##New Viewer Commands

IN ANY MODE

‘a’ - Switch to Animation mode

IN ANIMATION MODE

’s’ - Start and stop the animation running in the mesh view

‘q’ - Begin ray-tracing the animation, for 120 frames

’t’ - Switch between animation types (Starts with the wave equation, then the curvature flow by itself, then the full surface tension simulation)

‘d’ - Increases the displacement of vertices near the currently selected vertex by an amount proportional to that vertex’s displacement (helps make more reasonable displacements)

##Results

Our main results are videos showing our simulation of the wave equation running on a few different meshes with different material properties and initial conditions. Unfortunately, we didn’t have time to render particularly high quality videos, but what we have shows the wave equation dynamics and the water material working correctly. These videos are in Chris (ckaffine)’s handing directory, as is all of our code.

BVH Data:
On the dragon (dragon.dae from assignment3) -
Threads / Time to build (seconds)
1 - 0.3464
2 - 0.1879
4 - 0.1581
8 - 0.1218
Naive Implementation: 1 thread - 2.4 seconds 

