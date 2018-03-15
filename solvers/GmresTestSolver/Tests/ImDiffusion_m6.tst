<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady CG implicit diffusion </description>
    <executable>DiffusionSolver</executable>
    <parameters> ImDiffusion_m6.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_m6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-08"> 0.0020082 </value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-08"> 0.00277168 </value>
        </metric>
    </metrics>
</test>
