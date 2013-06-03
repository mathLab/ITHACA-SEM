<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LDG explicit diffusion, deformed quads, P=3</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LDG_deformed.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_2D_LDG_deformed.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00142525</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00922016</value>
        </metric>
    </metrics>
</test>




