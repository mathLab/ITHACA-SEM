<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>AcousticSolver</executable>
    <parameters>LEE_2DPulseAdv_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">LEE_2DPulseAdv_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p"   tolerance="1e-4">  6.77196</value>
            <value variable="rho" tolerance="1e-12"> 1.04459e-06</value>
            <value variable="rhou"  tolerance="1e-7">0.00191276</value>
            <value variable="rhov"  tolerance="1e-7">0.00195208</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p"   tolerance="1e-4">  30.1611</value>
            <value variable="rho" tolerance="1e-12"> 7.98548e-06</value>
            <value variable="rhou"  tolerance="1e-7">0.00713893</value>
            <value variable="rhov"  tolerance="1e-7">0.0098122</value>
        </metric>
    </metrics>
</test>
