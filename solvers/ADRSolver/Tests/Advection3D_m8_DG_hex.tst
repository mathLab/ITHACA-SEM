<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, hexahedra, order 4, P=12</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_m8_DG_hex.xml</parameters>
    <files>
        <file description="Session File">Advection3D_m8_DG_hex.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.08584e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.52297e-05</value>
        </metric>
    </metrics>
</test>

