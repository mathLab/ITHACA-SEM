<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> NekMesh with Periodic Boundary condition and Boundary Layer </description>
    <executable>NekMesh</executable>
    <parameters> -m peralign:dir=y:surf1=3:surf2=5 -m bl:surf=4,6:layers=4:r=3:nq=7 -m jac:list peralign_bl_cube.dat peralign_bl_cube-out.xml:xml:test </parameters>
    <files>
        <file description="Input File">peralign_bl_cube.dat</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
