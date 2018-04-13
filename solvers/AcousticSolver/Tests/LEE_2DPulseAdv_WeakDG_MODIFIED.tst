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
            <value variable="p"   tolerance="1e-4">6.77178</value>
            <value variable="rho" tolerance="1e-12">1.04863e-06</value>
            <value variable="ru"  tolerance="1e-7">0.0019274</value>
            <value variable="rv"  tolerance="1e-7">0.001952</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p"   tolerance="1e-4">30.1506</value>
            <value variable="rho" tolerance="1e-12">8.06966e-06</value>
            <value variable="ru"  tolerance="1e-7">0.00733383</value>
            <value variable="rv"  tolerance="1e-7">0.0098123</value>
        </metric>
    </metrics>
</test>
