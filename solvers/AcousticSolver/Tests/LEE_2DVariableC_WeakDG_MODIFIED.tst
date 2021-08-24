<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>AcousticSolver</executable>
    <parameters>LEE_2DVariableC_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">LEE_2DVariableC_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p"   tolerance="1e-4">43.2582</value>
            <value variable="rho" tolerance="1e-4">0.000580605</value>
            <value variable="rhou"  tolerance="1e-7">0.272808</value>
            <value variable="rhov"  tolerance="1e-7">0.0744645</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p"   tolerance="1e-4">53.3125</value>
            <value variable="rho" tolerance="1e-4">0.000835853</value>
            <value variable="rhou"  tolerance="1e-7">0.331939</value>
            <value variable="rhov"  tolerance="1e-7">0.120604</value>
        </metric>
    </metrics>
</test>
