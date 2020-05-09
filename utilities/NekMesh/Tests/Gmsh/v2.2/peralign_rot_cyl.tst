<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) tri mesh with rotated periodic boundary using peralign</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list -m peralign:rot="-PI/3.0":dir=z:surf1=5:surf2=4 peralign_rot_cyl.msh peralign_rot_cyl-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">peralign_rot_cyl.msh</file>
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
