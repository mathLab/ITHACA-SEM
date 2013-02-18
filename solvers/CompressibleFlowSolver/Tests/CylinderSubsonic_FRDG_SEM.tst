<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder, Dirichlet bcs, FRDG, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_FRDG_SEM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_FRDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">10.7745</value>
            <value variable="rhou" tolerance="1e-12">1.07744</value>
            <value variable="rhov" tolerance="1e-8">0.000452004</value>
            <value variable="E" tolerance="1e-12">0.131822</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.28742</value>
            <value variable="rhou" tolerance="1e-12">0.124093</value>
            <value variable="rhov" tolerance="1e-8">0.00381254</value>
            <value variable="E" tolerance="1e-12">0.0159379</value>
        </metric>
    </metrics>
</test>


