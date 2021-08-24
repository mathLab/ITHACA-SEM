<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady LDG explicit diffusion, periodic in x</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_2D_LDG_hybrid_m3_periodic.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_2D_LDG_hybrid_m3_periodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00170534</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00472021</value>
        </metric>
    </metrics>
</test>