<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D Helmholtz/Steady Diffusion Reaction P=7 </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz2D_modal.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_modal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.00888037 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.0101781 </value>
        </metric>
    </metrics>
</test>
