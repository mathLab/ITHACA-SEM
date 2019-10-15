<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.1) tet cube, convert to order 9</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list cube_tet_lin.msh cube_tet_lin-out.xml:xml:order=9:test</parameters>
    <files>
        <file description="Input File">cube_tet_lin.msh</file>
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
