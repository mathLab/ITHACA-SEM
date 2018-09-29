<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.0) linear prism with order 4 output</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list CubePrism_4_0.msh CubePrism_4_0.xml:xml:test:order=4</parameters>
    <files>
        <file description="Input File">CubePrism_4_0.msh</file>
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
