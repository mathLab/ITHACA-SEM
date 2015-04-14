<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, FRDG advection and LFRDG diffusion, GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_FRDG_LFRDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_FRDG_LFRDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">5.50779</value>
            <value variable="rhou" tolerance="1e-12">2016.32</value>
            <value variable="rhov" tolerance="1e-8">24.8084</value>
            <value variable="E" tolerance="1e-12">6.27023e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.33243</value>
            <value variable="rhou" tolerance="1e-12">87.3714</value>
            <value variable="rhov" tolerance="1e-8">21.6213</value>
            <value variable="E" tolerance="1e-12">287936</value>
        </metric>
    </metrics>
</test>


