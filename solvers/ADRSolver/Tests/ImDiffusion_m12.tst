<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion, order 3, P=12 </description>
    <executable>ADRSolver</executable>
    <parameters> ImDiffusion_m12.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_m12.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08"> 1.08903e-07 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08"> 1.48052e-06 </value>
        </metric>
    </metrics>
</test>
