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
            <value variable="p"   tolerance="1e-4">43.259</value>
            <value variable="rho" tolerance="1e-4">0.000580608</value>
            <value variable="ru"  tolerance="1e-7">0.272809</value>
            <value variable="rv"  tolerance="1e-7">0.0744649</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p"   tolerance="1e-4">53.3125</value>
            <value variable="rho" tolerance="1e-4">0.000835853</value>
            <value variable="ru"  tolerance="1e-7">0.331939</value>
            <value variable="rv"  tolerance="1e-7">0.120604</value>
        </metric>
    </metrics>
</test>
