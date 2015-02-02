<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh high-order prism cube, order 2</description>
    <executable>MeshConvert</executable>
    <parameters>-m jac:list CubePrism.msh CubePrism.xml:xml:test</parameters>
    <files>
        <file description="Input File">CubePrism.msh</file>
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
