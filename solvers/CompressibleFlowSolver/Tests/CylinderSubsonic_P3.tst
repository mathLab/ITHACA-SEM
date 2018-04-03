<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder P=3</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_P3.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_P3.xml</file>
        <file description="Restart File">CylinderSubsonic_P3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">42.8142</value>
            <value variable="rhou" tolerance="1e-12">4.28266</value>
            <value variable="rhov" tolerance="1e-12">0.0825261</value>
            <value variable="E" tolerance="1e-12">5.23742</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.32134</value>
            <value variable="rhou" tolerance="1e-12">0.199461</value>
            <value variable="rhov" tolerance="1e-12">0.114291</value>
            <value variable="E" tolerance="1e-12">0.159844</value>
        </metric>
    </metrics>
</test>


