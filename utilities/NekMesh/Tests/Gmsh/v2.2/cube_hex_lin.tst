<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) linear hex with order 7 output</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list cube_hex_lin.msh cube_hex_lin-out.xml:xml:test:order=7</parameters>
    <files>
        <file description="Input File">cube_hex_lin.msh</file>
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
