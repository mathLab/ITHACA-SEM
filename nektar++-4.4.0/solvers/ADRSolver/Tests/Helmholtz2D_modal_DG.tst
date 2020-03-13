<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D Helmholtz/Steady Diffusion Reaction (DG) P=7 </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz2D_modal_DG.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_modal_DG.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.011251 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.0355625 </value>
        </metric>
    </metrics>
</test>
