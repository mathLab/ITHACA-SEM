<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>MMF SWE solver, DG, P=4</description>
    <executable>ShallowWaterSolver</executable>
    <parameters>MMFSWERossbySphere.xml</parameters>
    <files>
        <file description="Session File">MMFSWERossbySphere.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="eta" tolerance="1e-5">0.000537129</value>
            <value variable="u" tolerance="1e-5">0.0157969</value>
            <value variable="v" tolerance="1e-5">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="eta" tolerance="1e-5"> 0.00524444</value>
            <value variable="u" tolerance="1e-5">0.0546249</value>
            <value variable="v" tolerance="1e-5">0.0</value>
        </metric>
    </metrics>
</test>


