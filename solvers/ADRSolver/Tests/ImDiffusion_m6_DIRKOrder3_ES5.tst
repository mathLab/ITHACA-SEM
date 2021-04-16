<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion </description>
    <executable>ADRSolver</executable>
    <parameters> ImDiffusion_m6_DIRKOrder3_ES5.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_m6_DIRKOrder3_ES5.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08"> 1.77579e-03 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08"> 2.05780e-02 </value>
        </metric>
    </metrics>
</test>
