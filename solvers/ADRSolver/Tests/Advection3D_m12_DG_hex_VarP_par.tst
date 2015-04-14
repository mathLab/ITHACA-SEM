<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, hexahedra, order 1, P=12,periodic bcs</description>
    <executable>ADRSolver</executable>
    <parameters>--use-metis Advection3D_m12_DG_hex_VarP.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Advection3D_m12_DG_hex_VarP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.99686e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-09">1.50105e-05</value>
        </metric>
    </metrics>
</test>
