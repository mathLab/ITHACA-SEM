<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, 2 hexahedra, covering all eDir1xxxDir2_Dir2xxxDir1 combinations</description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_DG_hex_faceRotation1221.xml Advection3D_DG_hex_faceRotation.xml</parameters>
    <files>
        <file description="Session File">Advection3D_DG_hex_faceRotation.xml</file>
        <file description="Mesh File">Advection3D_DG_hex_faceRotation1221.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00339959</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00174967</value>
        </metric>
    </metrics>
</test>
