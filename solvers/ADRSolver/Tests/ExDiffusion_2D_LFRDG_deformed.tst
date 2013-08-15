<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LFRDG explicit diffusion, deformed quads, P=3</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LFRDG_deformed.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_2D_LFRDG_deformed.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00165715</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0110807</value>
        </metric>
    </metrics>
</test>




