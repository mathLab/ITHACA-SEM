<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler 2D shocktube with mixed mesh and Laplacian AV</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_Lap.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_Lap.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000367418</value>
            <value variable="rhou" tolerance="1e-10">0.110602</value>
            <value variable="rhov" tolerance="1e-12">0.00309297</value>
            <value variable="E" tolerance="1e-6">999.27</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-11">0.0113659</value>
            <value variable="rhou" tolerance="1e-9">2.70304</value>
            <value variable="rhov" tolerance="1e-10">0.406501</value>
            <value variable="E" tolerance="1e-6">30654.9</value>
        </metric>
    </metrics>
</test>