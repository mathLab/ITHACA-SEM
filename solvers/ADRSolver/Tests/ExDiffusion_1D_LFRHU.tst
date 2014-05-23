<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>1D unsteady LFRHU explicit diffusion, order 3, P=2</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_1D_LFRHU.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_1D_LFRHU.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000999975</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00107628</value>
        </metric>
    </metrics>
</test>




