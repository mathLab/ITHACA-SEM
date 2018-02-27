<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, hexahedra, order 1, P=12,periodic bcs</description>
    <executable>ADRSolver</executable>
    <parameters>--use-scotch Advection3D_m12_DG_hex_periodic.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Advection3D_m12_DG_hex_periodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000497955</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-03">0.000612751</value>
        </metric>
    </metrics>
</test>
