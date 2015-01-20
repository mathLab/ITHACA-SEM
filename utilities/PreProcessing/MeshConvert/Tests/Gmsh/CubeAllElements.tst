<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh linear mesh of cube, all elements</description>
    <executable>MeshConvert</executable>
    <parameters>-m jac:list CubeAllElements.msh CubeAllElements.xml:xml:test</parameters>
    <files>
        <file description="Input File">CubeAllElements.msh</file>
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
