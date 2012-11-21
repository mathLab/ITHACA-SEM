<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion </description>
    <executable>ADRSolver</executable>
    <parameters> ImDiffusion_m6.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_m6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.00424554 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.032055 </value>
        </metric>
    </metrics>
</test>
