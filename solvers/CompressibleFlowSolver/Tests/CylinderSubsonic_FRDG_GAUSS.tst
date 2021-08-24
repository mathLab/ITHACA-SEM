<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder, Dirichlet bcs, FRDG, GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_FRDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_FRDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">10.7744</value>
            <value variable="rhou" tolerance="1e-12">1.07718</value>
            <value variable="rhov" tolerance="1e-8">0.00891406</value>
            <value variable="E" tolerance="1e-12">2.228e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.22526</value>
            <value variable="rhou" tolerance="1e-12">0.132264</value>
            <value variable="rhov" tolerance="1e-8">0.0468928</value>
            <value variable="E" tolerance="1e-12">253388</value>
        </metric>
    </metrics>
</test>


