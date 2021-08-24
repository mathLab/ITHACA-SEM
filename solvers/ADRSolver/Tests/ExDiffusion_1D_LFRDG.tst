<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady LFRDG explicit diffusion, order 3, P=2</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_1D_LFRDG.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_1D_LFRDG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00100248</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00104541</value>
        </metric>
    </metrics>
</test>




