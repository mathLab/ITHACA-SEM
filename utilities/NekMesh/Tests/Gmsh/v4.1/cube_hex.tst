<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.1) high-order hex cube</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list cube_hex.msh cube_hex-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">cube_hex.msh</file>
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
