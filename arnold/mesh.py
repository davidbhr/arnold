import gmsh
import numpy as np
from time import sleep
import os



def cylindrical_inclusion(mesh_file, d_cyl=30, l_cyl=300, r_outer=2000, length_factor=0.2):
    """
    Creates a spherical bulk mesh with a centered cylindrical inclusion.
    
    Args:
        mesh_file(str): File path to save mesh file (add .msh ending)
        d_cyl(float): Diameter of the spherical inclusion (in µm)
        l_cyl(float): Length of the spherical inclusion (in µm)
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
    """
    
    # create subfolders to mesh path if they do not exist --------------------------------------------------------------
    os.makedirs(os.path.dirname(mesh_file), exist_ok=True)
            
    # Build model (OpenCascade (for boolean operators) or built-in)
    model = gmsh.model
    factory = model.occ  # occ opencascade for boolean operations,  else (factory = model.geo)
    gmsh.initialize('', False)  # do not read in config files,  else (gmsh.initialize(sys.argv))
    gmsh.option.setNumber("General.Terminal", 1)
    model.add("TwoSpheres")  # Model Name
    
    # Mesh algorithm ---------------------------------------------------------------------------------------------------
    # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal, 7: BAMG, 8: DelQuad)
    gmsh.option.setNumber("Mesh.Algorithm", 2)
    
    #3D mesh algorithm (1: Delaunay, 4: Frontal, 5: Frontal Delaunay, 6: Frontal Hex, 7: MMG3D, 9: R-tree, 10: HXT)
    gmsh.option.setNumber("Mesh.Algorithm3D", 4)
    
    # Mesh size --------------------------------------------------------------------------------------------------------
    # Factor applied to all mesh element sizes, Default value: 1  
    gmsh.option.setNumber('Mesh.CharacteristicLengthFactor', length_factor)   
    
    # Element Size by Curvature (0=off 1= on, Include curvature into element mesh size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)
    
    # Extend computation of mesh element sizes from the boundaries into the interior 
    # (for 3D Delaunay, use 1: longest or 2: shortest surface edge length; 1 by default)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
    
     # Turn off 0 for dimensionless mesh, or if 1 set mesh size (so far only possible for points..)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    
    # Optimize the Mesh-------------------------------------------------------------------------------------------------
    # Optimize the mesh to improve the quality of tetrahedral elements, Default value: 1
    gmsh.option.setNumber("Mesh.Optimize", 1)

    # Optimize tetrahedra that have a quality below ...,  Default value: 0.3       
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.3)
    
    # Optimize the mesh using Netgen to improve the quality of tetrahedral elements, Default value: 0
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    
    
    # Mesh Output Format -----------------------------------------------------------------------------------------------
    # default=10  (1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos stat, 27: stl, 28: p3d, 30: mesh, 
    #              31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39: inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 
    #              49: neu, 50: matlab)
    gmsh.option.setNumber("Mesh.Format", 1)
    
    # Currently select version 2 (2.2) to apply appropriate boundary conditions
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    
    
    # Create Geometry --------------------------------------------------------------------------------------------------
    factory.addCylinder(-l_cyl/2*1e-6,0,0, l_cyl*1e-6,0,0, (d_cyl/2)*1e-6, tag=1, angle=2*np.pi) # cylindric inclusion
    factory.addSphere(0,0,0, r_outer*1e-6, tag=2, angle1=-np.pi/2, angle2=np.pi/2, angle3=2*np.pi) # bulk
        
    # Boolean cut outer Sphere
    # cut: boolean difference   
    # Input: objectDimTags, toolDimTags, tag, removeObject, removeTool      
    # Output: outDimTags, outDimTagsMap 
    factory.cut([(3,2)], [(3,1)], tag=3, removeObject=True, removeTool=True) # boolean cut
    factory.synchronize()
    
    # Define Physical Group (to apply boundary conditions later-on)
    model.addPhysicalGroup(3,[3], 4)  # (dim, tags, newtag)
    
    #Synchronize and generate mesh
    factory.synchronize()
    model.mesh.generate(3)
    
    # Save the mesh model
    gmsh.write(mesh_file)
    sleep(5)  # wait for file to be written to disk






def spherical_inclusion(mesh_file, r_inner=100, r_outer=20000, length_factor=0.2):
    """
    Creates a spherical bulk mesh with a centered spherical inclusion.
    
    Args:
        mesh_file(str): File path to save mesh file (add .msh ending)
        r_inner(float): Radius of the spherical inclusion (in µm)
        r_outer(float): Outer radius of the bulk mesh (in µm)
        length_factor(float): Mesh element size is determined by curvature and then multipled with this factor
    """
    # Build model (OpenCascade (for boolean operators) or built-in)
    model = gmsh.model
    factory = model.occ  # occ opencascade for boolean operations,  else (factory = model.geo)
    gmsh.initialize('', False)  # do not read in config files,  else (gmsh.initialize(sys.argv))
    gmsh.option.setNumber("General.Terminal", 1)
    model.add("TwoSpheres")  # Model Name
    
    # Mesh algorithm ---------------------------------------------------------------------------------------------------
    # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal, 7: BAMG, 8: DelQuad)
    gmsh.option.setNumber("Mesh.Algorithm", 2)
    
    #3D mesh algorithm (1: Delaunay, 4: Frontal, 5: Frontal Delaunay, 6: Frontal Hex, 7: MMG3D, 9: R-tree, 10: HXT)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    
    # Mesh size --------------------------------------------------------------------------------------------------------
    # Factor applied to all mesh element sizes, Default value: 1  
    gmsh.option.setNumber('Mesh.CharacteristicLengthFactor', length_factor)   
    
    # Element Size by Curvature (0=off 1= on, Include curvature into element mesh size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)
    
    # Extend computation of mesh element sizes from the boundaries into the interior 
    # (for 3D Delaunay, use 1: longest or 2: shortest surface edge length; 1 by default)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
    
     # Turn off 0 for dimensionless mesh, or if 1 set mesh size (so far only possible for points..)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    
    # Optimize the Mesh-------------------------------------------------------------------------------------------------
    # Optimize the mesh to improve the quality of tetrahedral elements, Default value: 1
    gmsh.option.setNumber("Mesh.Optimize", 1)

    # Optimize tetrahedra that have a quality below ...,  Default value: 0.3       
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.3)
    
    # Optimize the mesh using Netgen to improve the quality of tetrahedral elements, Default value: 0
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    
    
    # Mesh Output Format -----------------------------------------------------------------------------------------------
    # default=10  (1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos stat, 27: stl, 28: p3d, 30: mesh, 
    #              31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39: inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 
    #              49: neu, 50: matlab)
    gmsh.option.setNumber("Mesh.Format", 1)
    
    # Currently select version 2 (2.2) to apply appropriate boundary conditions
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    
    # Create Geometry---------------------------------------------------------------------------------------------------
    factory.addSphere(0,0,0, r_inner*1e-6, tag=1, angle1=-np.pi/2, angle2=np.pi/2, angle3=2*np.pi)  # inclusion
    factory.addSphere(0,0,0, r_outer*1e-6, tag=2, angle1=-np.pi/2, angle2=np.pi/2, angle3=2*np.pi)  # bulk
    
    # Boolean cut outer Sphere
    # cut: boolean difference   
    # Input: objectDimTags, toolDimTags, tag, removeObject, removeTool      
    # Output: outDimTags, outDimTagsMap 
    factory.cut([(3,2)], [(3,1)], tag=3, removeObject=True, removeTool=True)  # boolean cut
    factory.synchronize()
    
    # Define Physical Group (to apply boundary conditions later-on)
    model.addPhysicalGroup(3,[3], 4)  # (dim, tags, newtag)
    
    #Synchronize and generate mesh
    factory.synchronize()
    model.mesh.generate(3)
    
    # Save the mesh model
    gmsh.write(mesh_file)
    sleep(1)  # wait for file to be written to disk
    
      

def show_mesh(mesh_file):
    """
    Opens and shows mesh file in Gmsh
    
    Args:
        mesh_file(str): File path to load mesh file
    """
    gmsh.initialize('', False)  # do not read in config files
    gmsh.open(mesh_file)
    gmsh.fltk.run()
