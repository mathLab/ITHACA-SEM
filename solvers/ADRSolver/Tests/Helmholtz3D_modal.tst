<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction P=7 With all types of BCs </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_modal.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_modal.xml</file>
    </files>
    <metrics>
      <metric type="L2" id="1">
            <value variable="w" tolerance="1e-12">  2.54892e-06 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="w" tolerance="1e-12"> 5.41863e-06 </value>
        </metric>
    </metrics>
</test>
