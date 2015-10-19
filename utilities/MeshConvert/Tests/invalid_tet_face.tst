<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Tet with an invalid face curvature</description>
    <executable>MeshConvert</executable>
    <parameters>-m jac:list invalid_tet_face.xml tet_Face.xml:xml:test</parameters>
    <files>
        <file description="Input File">invalid_tet_face.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
